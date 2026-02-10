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
from molvs.fragment import LargestFragmentChooser
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

        # Step 3: Select largest fragment (removes counterions/salts)
        lfc = LargestFragmentChooser(prefer_organic=True)
        mol_parent = lfc.choose(mol_std)

        frags_after = len(Chem.GetMolFrags(mol_parent))
        result["num_fragments_removed"] = frags_before - frags_after

        if result["num_fragments_removed"] > 0:
            result["standardization_notes"].append(
                f"Removed {result['num_fragments_removed']} salt/fragment(s)"
            )

        # Step 4: Neutralize charges
        original_charge = Chem.GetFormalCharge(mol_parent)
        uc = Uncharger()
        mol_neutral = uc.uncharge(mol_parent)
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
