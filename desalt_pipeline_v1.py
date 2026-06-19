"""
Desalting & Neutralization Pipeline for Mass Spectrometry Databases
===================================================================
Reads a table with SMILES + chemical metadata, produces desalted/neutralized
parent structures with updated molecular properties.

Requirements:
    pip install rdkit molvs pandas openpyxl

Usage:
    python desalt_pipeline.py input.csv output.csv --smiles-col SMILES --delimiter ","
    python desalt_pipeline.py input.xlsx output.xlsx --smiles-col SMILES
"""

import argparse
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, inchi, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from molvs import Standardizer
from molvs.charge import Uncharger
import logging
import os
import warnings

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

EXCEL_EXTENSIONS = (".xlsx", ".xls")

TRANSITION_METALS = {
    "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
    "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
    "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
    "La","Ce","Pr","Nd","Sm","Eu","Gd","Tb","Dy",
}


# ── I/O helpers ────────────────────────────────────────────────────────────

def detect_encoding(filepath):
    """Try UTF-8 first; fall back to latin-1 (which accepts any byte sequence)."""
    with open(filepath, "rb") as f:
        raw = f.read()
    try:
        raw.decode("utf-8")
        return "utf-8"
    except UnicodeDecodeError:
        return "latin-1"


def find_header_row(filepath, target_col):
    """Scan an Excel file for the row containing target_col and return its index."""
    df_raw = pd.read_excel(filepath, header=None, dtype=str)
    for i, row in df_raw.iterrows():
        if target_col in row.values:
            return i
    return None


# ── Core standardization functions ──────────────────────────────────────────

def _choose_largest_fragment(mol):
    """Select the parent fragment from a multi-fragment (salt) molecule.

    Ranks fragments by HEAVY-ATOM count, then prefers a formally positive fragment
    (keeps quaternary-ammonium actives on a heavy-atom tie, e.g. DADMAC), then exact
    MW; organic fragments are preferred over inorganic. Returns
    (parent_mol, ambiguous), where `ambiguous` is True when the top two fragments are
    within one heavy atom — a near-tie worth flagging for manual review.

    Replaces MolVS/RDKit `LargestFragmentChooser`, whose 'largest' metric counts
    implicit hydrogens — so an H-rich aliphatic counterion (e.g. the diisopropylamine
    of a 2,4-D salt: 22 atoms-with-H) wrongly outranks an H-poor halogenated active
    (2,4-D: 19 atoms-with-H but 13 heavy atoms vs 7). Heavy-atom ranking keeps the
    active; legitimate desalts (Na/K/Cl/sulfate counterions, etc.) are unaffected.
    """
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) <= 1:
        return mol, False
    organic = [f for f in frags if any(a.GetAtomicNum() == 6 for a in f.GetAtoms())]
    pool = organic or list(frags)
    ranked = sorted(
        pool,
        key=lambda f: (f.GetNumHeavyAtoms(),
                       1 if Chem.GetFormalCharge(f) > 0 else 0,
                       Descriptors.ExactMolWt(f)),
        reverse=True,
    )
    top = ranked[0]
    # Near-tie flag: a *different* fragment within one heavy atom of the chosen one
    # (ignore identical duplicate fragments, e.g. a di-salt of the same acid).
    top_smi = Chem.MolToSmiles(top)
    rivals = [f for f in ranked[1:] if Chem.MolToSmiles(f) != top_smi]
    ambiguous = bool(rivals) and (top.GetNumHeavyAtoms() - rivals[0].GetNumHeavyAtoms() <= 1)
    return top, ambiguous


def balance_permanent_cation(mol):
    """Restore zwitterions that Uncharger leaves at net +1.

    MolVS `Uncharger` neutralizes by adding/removing protons, but a PERMANENT
    cation (quaternary N+, sulfonium S+, aromatic n+) has no proton to remove. For
    a zwitterion that has such a cation AND a free carboxylic acid, Uncharger
    protonates the carboxylate (COO- -> COOH) and the molecule comes out at net +1
    instead of the correct neutral zwitterion — giving a one-proton-too-high exact
    mass, a charged formula, and a protonated InChIKey (…-O instead of …-N). This
    breaks adduct m/z and database matching for carnitine, acylcarnitines, betaines,
    etc.

    If `mol` is net-positive due to a permanent cation AND has free carboxylic
    acid group(s), deprotonate one carboxylate per +1 charge to restore the neutral
    zwitterion. Returns (mol, changed). True cations with no free -COOH
    (choline, acetylcholine — an ester, thiamine) and neutrals are returned
    unchanged.
    """
    if Chem.GetFormalCharge(mol) <= 0:
        return mol, False
    has_permanent_cation = any(
        atom.GetFormalCharge() > 0 and (
            # quaternary N+ = 4 HEAVY (non-H) substituents. GetDegree() excludes
            # implicit/explicit H, so this excludes ammonium (NH4+, degree 0) and
            # protonated amines (R-NH3+/R2NH2+/R3NH+, degree <4) — both counterions
            # that Uncharger neutralizes, not permanent cations.
            (atom.GetAtomicNum() == 7 and atom.GetDegree() == 4) or         # quaternary N+
            (atom.GetAtomicNum() == 16 and atom.GetFormalCharge() > 0) or   # sulfonium S+
            (atom.GetAtomicNum() == 7 and atom.GetIsAromatic()             # aromatic n+
             and atom.GetFormalCharge() > 0)
        )
        for atom in mol.GetAtoms()
    )
    if not has_permanent_cation:
        return mol, False
    matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3](=O)[OX2H1]"))
    if not matches:                       # no free carboxylate to balance it -> leave as-is
        return mol, False

    rw = Chem.RWMol(mol)
    remaining = Chem.GetFormalCharge(mol)
    for match in matches:
        if remaining <= 0:
            break
        oh = rw.GetAtomWithIdx(match[2])  # the carboxyl -OH oxygen
        n_h = oh.GetTotalNumHs()
        if n_h > 0:
            oh.SetNoImplicit(True)
            oh.SetNumExplicitHs(n_h - 1)
            oh.SetFormalCharge(-1)
            remaining -= 1
    new_mol = rw.GetMol()
    try:
        Chem.SanitizeMol(new_mol)
    except Exception:
        return mol, False                 # never emit a broken molecule
    return new_mol, True


def standardize_mol(smiles: str) -> dict:
    """
    Full standardization pipeline for a single SMILES string.
    
    Steps:
        1. Parse SMILES
        2. Standardize (normalize functional groups, canonical tautomers)
        3. Remove salts / select largest fragment
        4. Neutralize charges
        5. Compute updated properties
    
    Returns dict with MS-ready SMILES, InChIKey, formula, exact mass, etc.
    """
    result = {
        "msready_smiles": None,
        "msready_inchikey": None,
        "msready_formula": None,
        "msready_exact_mass": None,
        "msready_monoisotopic_mass": None,
        "num_fragments_removed": 0,
        "charge_changed": False,
        "standardization_notes": [],
    }

    # Step 1: Parse
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        result["standardization_notes"] = "FAILED: Could not parse SMILES"
        return result

    # Check for transition metals — skip desalting for organometallics
    atom_symbols = {atom.GetSymbol() for atom in mol.GetAtoms()}
    found_metals = atom_symbols & TRANSITION_METALS
    if found_metals:
        result["msready_smiles"] = Chem.MolToSmiles(mol, canonical=True)
        result["msready_inchikey"] = inchi.MolToInchiKey(mol)
        result["msready_formula"] = rdMolDescriptors.CalcMolFormula(mol)
        result["msready_exact_mass"] = round(Descriptors.ExactMolWt(mol), 6)
        result["msready_monoisotopic_mass"] = round(Descriptors.ExactMolWt(mol), 6)
        result["standardization_notes"] = f"SKIPPED: Contains transition metal ({', '.join(sorted(found_metals))}); manual review recommended"
        return result

    try:
        # Step 2: MolVS Standardizer (normalizes nitro groups, charges, etc.)
        s = Standardizer()
        mol_std = s.standardize(mol)

        # Count fragments before salt removal
        frags_before = len(Chem.GetMolFrags(mol_std))

        # Step 3: Select the parent fragment (removes counterions/salts) by
        # HEAVY-ATOM count — not LargestFragmentChooser, whose implicit-H counting
        # mis-picks H-rich counterions of H-poor actives (e.g. 2,4-D amine salts).
        mol_parent, desalt_ambiguous = _choose_largest_fragment(mol_std)

        frags_after = len(Chem.GetMolFrags(mol_parent))
        result["num_fragments_removed"] = frags_before - frags_after

        if result["num_fragments_removed"] > 0:
            result["standardization_notes"].append(
                f"Removed {result['num_fragments_removed']} salt/fragment(s)"
            )
        if desalt_ambiguous:
            result["standardization_notes"].append(
                "Desalt ambiguous (near-tie fragments) — review"
            )

        # Step 4: Neutralize charges
        original_charge = Chem.GetFormalCharge(mol_parent)
        uc = Uncharger()
        mol_neutral = uc.uncharge(mol_parent)

        # Step 4b: Re-balance zwitterions Uncharger left at net +1 (permanent
        # cation + free carboxylate) so the MS-Ready structure is truly neutral.
        mol_neutral, zwitterion_fixed = balance_permanent_cation(mol_neutral)
        if zwitterion_fixed:
            result["standardization_notes"].append(
                "Re-balanced permanent-cation zwitterion (restored carboxylate)"
            )

        final_charge = Chem.GetFormalCharge(mol_neutral)

        if original_charge != final_charge:
            result["charge_changed"] = True
            result["standardization_notes"].append(
                f"Neutralized charge: {original_charge:+d} → {final_charge:+d}"
            )

        # Step 5: Compute properties on MS-ready structure
        result["msready_smiles"] = Chem.MolToSmiles(mol_neutral, canonical=True)
        result["msready_inchikey"] = inchi.MolToInchiKey(mol_neutral)
        result["msready_formula"] = rdMolDescriptors.CalcMolFormula(mol_neutral)
        result["msready_exact_mass"] = round(Descriptors.ExactMolWt(mol_neutral), 6)
        result["msready_monoisotopic_mass"] = round(
            Descriptors.ExactMolWt(mol_neutral), 6
        )

        if not result["standardization_notes"]:
            result["standardization_notes"].append("No changes needed")

    except Exception as e:
        result["standardization_notes"].append(f"ERROR: {str(e)}")

    # Collapse notes list to string
    result["standardization_notes"] = "; ".join(result["standardization_notes"])
    return result


# ── Batch processing ────────────────────────────────────────────────────────

def process_dataframe(
    df: pd.DataFrame,
    smiles_col: str = "SMILES",
) -> pd.DataFrame:
    """
    Process an entire DataFrame of chemical records.
    
    Args:
        df: Input DataFrame with at minimum a SMILES column
        smiles_col: Name of the column containing SMILES strings
    
    Returns:
        DataFrame with original columns + new MS-ready structure columns
    """
    if smiles_col not in df.columns:
        raise ValueError(
            f"Column '{smiles_col}' not found. Available: {list(df.columns)}"
        )

    logger.info(f"Processing {len(df)} compounds...")

    # Track statistics
    n_total = len(df)
    n_failed = 0
    n_desalted = 0
    n_neutralized = 0

    results = []
    for idx, row in df.iterrows():
        smi = row[smiles_col]

        if pd.isna(smi) or str(smi).strip() == "":
            results.append(
                {
                    "msready_smiles": None,
                    "msready_inchikey": None,
                    "msready_formula": None,
                    "msready_exact_mass": None,
                    "msready_monoisotopic_mass": None,
                    "num_fragments_removed": 0,
                    "charge_changed": False,
                    "standardization_notes": "FAILED: Empty SMILES",
                }
            )
            n_failed += 1
            continue

        res = standardize_mol(str(smi).strip())
        results.append(res)

        if res["msready_smiles"] is None:
            n_failed += 1
        else:
            if res["num_fragments_removed"] > 0:
                n_desalted += 1
            if res["charge_changed"]:
                n_neutralized += 1

        # Progress logging
        if (idx + 1) % 1000 == 0:
            logger.info(f"  Processed {idx + 1}/{n_total}...")

    # Merge results back
    results_df = pd.DataFrame(results)
    output = pd.concat([df.reset_index(drop=True), results_df], axis=1)

    # Summary
    logger.info("=" * 60)
    logger.info(f"SUMMARY")
    logger.info(f"  Total compounds:    {n_total}")
    logger.info(f"  Successfully parsed: {n_total - n_failed}")
    logger.info(f"  Failed to parse:     {n_failed}")
    logger.info(f"  Desalted:            {n_desalted}")
    logger.info(f"  Neutralized:         {n_neutralized}")
    logger.info("=" * 60)

    return output


# ── CLI entry point ─────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Desalt and neutralize chemical structures in a data table."
    )
    parser.add_argument("input_file", help="Input CSV/TSV/Excel file path")
    parser.add_argument("output_file", help="Output CSV/TSV/Excel file path")
    parser.add_argument(
        "--smiles-col",
        default="SMILES",
        help="Name of the SMILES column (default: SMILES)",
    )
    parser.add_argument(
        "--delimiter",
        default=",",
        help="File delimiter: ',' for CSV, '\\t' for TSV (default: ',')",
    )
    parser.add_argument(
        "--encoding",
        default="auto",
        help="File encoding for text files: 'auto' to detect, or explicit e.g. 'latin-1' (default: auto)",
    )
    parser.add_argument(
        "--header-row",
        type=int,
        default=None,
        help="Row number (0-indexed) containing column headers in Excel files (default: auto-detect)",
    )

    args = parser.parse_args()

    # Handle tab delimiter from CLI
    delim = "\t" if args.delimiter in ("\\t", "tab", "TAB") else args.delimiter

    # Read input
    logger.info(f"Reading {args.input_file}...")
    ext_in = os.path.splitext(args.input_file)[1].lower()
    if ext_in in EXCEL_EXTENSIONS:
        header_row = args.header_row
        if header_row is None:
            header_row = find_header_row(args.input_file, args.smiles_col)
            if header_row is not None:
                logger.info(f"  Auto-detected header at row {header_row}")
            else:
                logger.warning(f"  Could not find '{args.smiles_col}' column; using row 0 as header")
                header_row = 0
        df = pd.read_excel(args.input_file, header=header_row, dtype=str)
    else:
        encoding = detect_encoding(args.input_file) if args.encoding == "auto" else args.encoding
        logger.info(f"  Using encoding: {encoding}")
        df = pd.read_csv(args.input_file, sep=delim, dtype=str, encoding=encoding)
    logger.info(f"  Found {len(df)} rows, columns: {list(df.columns)}")

    # Process
    output = process_dataframe(df, smiles_col=args.smiles_col)

    # Write output
    ext_out = os.path.splitext(args.output_file)[1].lower()
    if ext_out in EXCEL_EXTENSIONS:
        output.to_excel(args.output_file, index=False)
    else:
        out_delim = "\t" if ext_out == ".tsv" else ","
        output.to_csv(args.output_file, sep=out_delim, index=False)
    logger.info(f"Saved results to {args.output_file}")


if __name__ == "__main__":
    main()
