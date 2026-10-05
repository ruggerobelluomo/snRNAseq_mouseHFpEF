# Private inputs (not included)

The scripts read these study-specific files, which are not part of the public release:

| File | Read by | Columns / format |
|---|---|---|
| `animal_notes.csv` | `notebooks/_common.R` | `kind` (`qc_flag` or `biological_note`), `animal`, `note` |
| `candidate_genes.txt` | `notebooks/05_EC_DE_analysis.R`, `05c_...R` | one gene symbol per line |
| `lab_findings.csv` | `notebooks/09_HFpEF_model_validation.R` | `figure`, `measure`, `reported` |
| `lab_concordance_map.csv` | `notebooks/09_HFpEF_model_validation.R` | `phenotype`, `lab_finding`, `basis`, `lab_direction`, `lab_contrasts`, `readout`, ... |
