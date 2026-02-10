# Desalting & Neutralization Pipeline for Mass Spectrometry Databases

A Python pipeline that reads chemical data tables containing SMILES structures and produces desalted, neutralized parent structures with updated molecular properties. Designed for curating and standardizing compound libraries used in LC-HRMS and GC-HRMS non-targeted analysis workflows.

## Background

Mass spectrometry reference databases often contain compounds registered as salt forms (e.g., metformin hydrochloride, ibuprofen sodium, PFOS potassium salt). When matching experimental features to database entries by exact mass or molecular formula, you need the **parent (free base/free acid)** form — not the salt. This pipeline automates that conversion across large compound tables.

The pipeline applies a four-step standardization process to each structure:

1. **Standardize** — Normalize functional group representations and canonical tautomers (e.g., consistent nitro group notation, charge separation patterns)
2. **Desalt** — Identify and remove counterions by selecting the largest organic fragment (e.g., strip Na⁺, K⁺, Cl⁻, sulfate)
3. **Neutralize** — Remove formal charges where chemically appropriate (e.g., carboxylate⁻ → carboxylic acid, ammonium⁺ → amine)
4. **Recalculate properties** — Generate canonical SMILES, InChIKey, molecular formula, and exact monoisotopic mass on the parent structure

## Requirements

- Python ≥ 3.9
- Dependencies:

```bash
pip install rdkit molvs pandas openpyxl
```

> **Note:** On some systems, RDKit is installed via conda: `conda install -c conda-forge rdkit`

## File Overview

| File | Description |
|---|---|
| `desalt_pipeline.py` | Main pipeline module (CLI tool + importable library) |
| `demo_desalt.py` | Demo script with 10 example compounds showing typical use cases |
| `desalted_output.csv` | Example output from the demo run |

## Usage

### Option 1: Command Line

Process a CSV, TSV, or Excel file directly from the terminal:

```bash
# Basic usage (assumes SMILES column is named "SMILES")
python desalt_pipeline.py input.csv output.csv

# Specify a different SMILES column name
python desalt_pipeline.py input.csv output.csv --smiles-col "canonical_smiles"

# Tab-delimited input
python desalt_pipeline.py input.tsv output.tsv --delimiter "\t"

# Excel input and output
python desalt_pipeline.py input.xlsx output.xlsx

# Excel input, CSV output
python desalt_pipeline.py input.xlsx output.csv

# Excel file with branding/header rows (auto-detected)
python desalt_pipeline.py input.xlsx output.xlsx

# Excel file with manually specified header row (0-indexed)
python desalt_pipeline.py input.xlsx output.xlsx --header-row 13

# Explicit encoding for non-UTF-8 text files
python desalt_pipeline.py input.txt output.tsv --delimiter "\t" --encoding latin-1
```

**CLI Arguments:**

| Argument | Required | Default | Description |
|---|---|---|---|
| `input_file` | Yes | — | Path to input CSV/TSV/Excel file |
| `output_file` | Yes | — | Path for output file (format detected from extension) |
| `--smiles-col` | No | `SMILES` | Column name containing SMILES strings |
| `--delimiter` | No | `,` | Field delimiter for text files (`,` for CSV, `\t` for TSV) |
| `--encoding` | No | `auto` | Text file encoding (`auto` tries UTF-8 then latin-1, or specify explicitly) |
| `--header-row` | No | auto-detect | Row number (0-indexed) containing column headers in Excel files; auto-detects by searching for `--smiles-col` |

### Option 2: Python Import

Use `process_dataframe()` directly in your own scripts or notebooks:

```python
import pandas as pd
from desalt_pipeline import process_dataframe

# Read your compound table
df = pd.read_csv("my_database.csv")

# Run the pipeline
result = process_dataframe(df, smiles_col="SMILES")

# Inspect results
print(result[["compound_name", "SMILES", "msready_smiles", "msready_formula", 
              "msready_exact_mass", "standardization_notes"]])

# Save
result.to_csv("my_database_desalted.csv", index=False)
```

### Option 3: Single Compound

Standardize a single SMILES string:

```python
from desalt_pipeline import standardize_mol

# Ibuprofen sodium salt
result = standardize_mol("[Na+].CC(C)Cc1ccc(cc1)C(C)C([O-])=O")

print(result)
# {
#   'msready_smiles': 'CC(C)Cc1ccc(C(C)C(=O)O)cc1',
#   'msready_inchikey': 'HEFNNWSXXWATRW-UHFFFAOYSA-N',
#   'msready_formula': 'C13H18O2',
#   'msready_exact_mass': 206.13068,
#   'msready_monoisotopic_mass': 206.13068,
#   'num_fragments_removed': 1,
#   'charge_changed': True,
#   'standardization_notes': 'Removed 1 salt/fragment(s); Neutralized charge: -1 → +0'
# }
```

### Option 4: Integration with R via reticulate

For integration with existing R/XCMS workflows:

```r
library(reticulate)

# Point to your Python environment
use_python("/path/to/python")  # or use_condaenv("myenv")

# Import the module
desalt <- import_from_path("desalt_pipeline", path = "/path/to/pipeline/")

# Load your data in R, convert, process, convert back
library(readr)
df_r <- read_csv("compounds.csv")
df_py <- r_to_py(df_r)

result_py <- desalt$process_dataframe(df_py, smiles_col = "SMILES")
result_r <- py_to_r(result_py)

write_csv(result_r, "compounds_desalted.csv")
```

## Output Columns

The pipeline appends these columns to your original data:

| Column | Type | Description |
|---|---|---|
| `msready_smiles` | str | Canonical SMILES of the desalted, neutralized MS-ready structure |
| `msready_inchikey` | str | InChIKey of the MS-ready structure (useful for deduplication) |
| `msready_formula` | str | Molecular formula of the MS-ready structure |
| `msready_exact_mass` | float | Exact monoisotopic mass of the MS-ready structure (Da) |
| `msready_monoisotopic_mass` | float | Monoisotopic mass (same as exact mass; included for compatibility) |
| `num_fragments_removed` | int | Number of salt/counterion fragments stripped |
| `charge_changed` | bool | Whether formal charge was neutralized |
| `standardization_notes` | str | Human-readable log of all transformations applied |

## Example Transformations

| Compound | Input SMILES | Parent SMILES | Notes |
|---|---|---|---|
| Metformin HCl | `CN(C)C(=N)NC(=N)N.[H]Cl` | `CN(C)C(=N)NC(=N)N` | Removed HCl counterion |
| Ibuprofen sodium | `[Na+].CC(C)Cc1ccc(cc1)C(C)C([O-])=O` | `CC(C)Cc1ccc(C(C)C(=O)O)cc1` | Removed Na⁺, neutralized carboxylate |
| PFOS potassium salt | `[K+].OS(=O)(=O)C(F)(F)...C(F)(F)F` | `O=S(=O)(O)C(F)(F)...C(F)(F)F` | Removed K⁺ counterion |
| Sodium glutamate | `[Na+].OC(=O)C(N)CCC([O-])=O` | `NC(CCC(=O)O)C(=O)O` | Removed Na⁺, neutralized carboxylate |
| Caffeine | `Cn1c(=O)c2c(ncn2C)n(C)c1=O` | `Cn1c(=O)c2c(ncn2C)n(C)c1=O` | No changes needed |
| Invalid SMILES | `NOT_A_SMILES` | `None` | Gracefully flagged as failed |

## How the Standardization Steps Work

### Step 1: MolVS Standardizer

Applies a series of transformations to normalize chemical representations:

- Disconnects metal-organic bonds
- Normalizes charge-separated functional groups (e.g., `[N+]([O-])=O` → consistent nitro notation)
- Applies canonical tautomer selection
- Standardizes hypervalent atoms

### Step 2: Largest Fragment Selection

Uses `LargestFragmentChooser` with `prefer_organic=True`:

- Splits the molecule at disconnected fragment boundaries (`.` in SMILES)
- Selects the largest organic fragment by heavy atom count
- Removes common counterions: Na⁺, K⁺, Ca²⁺, Cl⁻, Br⁻, sulfate, phosphate, etc.
- Prefers carbon-containing fragments over inorganic ones

### Step 3: Charge Neutralization

Uses `Uncharger` to protonate/deprotonate to the neutral form:

- Carboxylates (R-COO⁻) → Carboxylic acids (R-COOH)
- Ammonium (R-NH₃⁺) → Amines (R-NH₂)
- Phenolates (Ar-O⁻) → Phenols (Ar-OH)
- Preserves permanently charged species (e.g., quaternary ammonium)

### Step 4: Property Calculation

All properties are computed on the final parent structure:

- **Canonical SMILES** — RDKit canonical form for consistent representation
- **InChIKey** — IUPAC International Chemical Identifier hash (27 characters); first 14 characters represent the connectivity layer, useful for deduplication across salt forms
- **Molecular formula** — Hill system notation
- **Exact monoisotopic mass** — Calculated from most abundant isotope of each element

## Error Handling

The pipeline is designed to process large databases robustly:

- **Invalid SMILES** — Flagged with `"FAILED: Could not parse SMILES"` in the notes column; other columns set to `None`
- **Empty/missing SMILES** — Flagged with `"FAILED: Empty SMILES"`
- **Unexpected errors** — Caught per-compound with error message in the notes column; does not halt batch processing
- **Progress logging** — Logs every 1,000 compounds during batch processing
- **Summary statistics** — Printed at the end showing total, successful, failed, desalted, and neutralized counts

## Performance Considerations

- Processing speed is approximately **500–1,000 compounds/second** on a modern CPU
- Memory usage scales linearly with the number of compounds
- For very large databases (>1M compounds), consider chunked processing:

```python
import pandas as pd
from desalt_pipeline import process_dataframe

chunks = pd.read_csv("large_database.csv", chunksize=50000, dtype=str)
results = []
for chunk in chunks:
    results.append(process_dataframe(chunk, smiles_col="SMILES"))

output = pd.concat(results, ignore_index=True)
output.to_csv("large_database_desalted.csv", index=False)
```

## Common Use Cases

### Building an MS Reference Library

```python
# Pull structures from a source database, desalt, then compute exact masses
# for matching against experimental m/z features
result = process_dataframe(df, smiles_col="SMILES")

# Use msready_exact_mass for [M+H]+ matching
result["mz_MpH"] = result["msready_exact_mass"] + 1.007276
result["mz_MmH"] = result["msready_exact_mass"] - 1.007276
result["mz_MpNa"] = result["msready_exact_mass"] + 22.989218
```

### Deduplicating Across Salt Forms

```python
# Multiple salt forms of the same drug collapse to one parent InChIKey
result = process_dataframe(df, smiles_col="SMILES")
deduplicated = result.drop_duplicates(subset="msready_inchikey", keep="first")
print(f"Reduced {len(result)} entries to {len(deduplicated)} unique parent structures")
```

### QC Check on an Existing Database

```python
# Identify which entries in your database are salt forms
result = process_dataframe(df, smiles_col="SMILES")
salts = result[result["num_fragments_removed"] > 0]
charged = result[result["charge_changed"] == True]
failed = result[result["msready_smiles"].isna()]

print(f"Salt forms found: {len(salts)}")
print(f"Charged species neutralized: {len(charged)}")
print(f"Failed to parse: {len(failed)}")
```

## Limitations and Caveats

- **Zwitterions**: Amino acids and other zwitterionic compounds may not neutralize to the expected form in all cases, as the "neutral" form is pH-dependent. The pipeline produces the fully protonated/deprotonated uncharged form.
- **Metallic complexes**: Organometallic compounds or metal chelates may not desalt correctly if the metal is part of the active structure rather than a counterion.
- **Stereochemistry**: The pipeline preserves stereochemistry (E/Z, R/S) from the input SMILES but does not validate or correct it.
- **Tautomers**: MolVS applies canonical tautomer selection, but the "correct" tautomer can be context-dependent. For sensitive applications, manual review of tautomer-rich compounds (e.g., keto-enol, amide-imidic acid) is recommended.
- **Multi-component salts**: For 2:1 or 1:2 salt complexes (e.g., calcium salts of diprotic acids), the pipeline selects only one copy of the largest fragment.

## Dependencies

| Package | Version | Purpose |
|---|---|---|
| [RDKit](https://www.rdkit.org/) | ≥ 2023.03 | Core cheminformatics: SMILES parsing, property calculation, InChIKey generation |
| [MolVS](https://molvs.readthedocs.io/) | ≥ 0.1.1 | Molecule Validation & Standardization: normalization, fragment selection, charge neutralization |
| [pandas](https://pandas.pydata.org/) | ≥ 1.5 | Data table I/O and manipulation |
| [openpyxl](https://openpyxl.readthedocs.io/) | ≥ 3.0 | Excel (.xlsx) file reading and writing (required only for Excel I/O) |

## License

MIT — free for academic and commercial use.
