# BRD4 Selective Inhibitor Discovery via Integrated Virtual Screening

A computational drug discovery pipeline combining QSAR pre-screening, molecular docking, and selectivity-guided candidate selection to identify BRD4-selective inhibitor candidates from approximately one million ZINC22 lead-like compounds.

This repository accompanies the PHM5013 Precision Drug Discovery and Pharmacogenomics group project (2026).

---

## Project Overview

The pipeline integrates three complementary approaches:

1. **Ligand-based screening** — Random Forest QSAR models trained on ChEMBL bioactivity data (BRD4 test R² = 0.779; BRD2 as qualitative selectivity filter)
2. **Structure-based docking** — Validated AutoDock Vina protocols for both BRD4 (PDB: 3MXF) and BRD2 (PDB: 7LAH), with retrospective ROC/EF enrichment analysis (BRD4 AUC = 0.864 with RF-Score-VS; BRD2 AUC = 0.942)
3. **Molecular dynamics** — 100 ns triplicate MD simulations with MM-GBSA binding free energy analysis

Three candidate compounds (LIG_1629, LIG_1610, LIG_0061) were identified with strong BRD4 docking affinities (−7.2 to −8.2 kcal/mol) and preliminary BRD4/BRD2 selectivity differentials of 1.8–2.4 kcal/mol.

---

## Repository Structure

```
.
├── 01_BRD4_QSAR_MODEL/             QSAR model development and validation (BRD4)
│   └── QSAR_BRD4.ipynb
│
├── 02_BRD2_QSAR_MODEL/             QSAR model development and validation (BRD2)
│   └── QSAR_BRD2.ipynb
│
├── 03_BRD4_VS_1M/                  HPC virtual screening of ~1M ZINC22 compounds
│   ├── hpc_mega_screen.py          Main screening script (RF model + AD assessment)
│   ├── run_vs_brd4.sh              PBS submission script
│   └── csv_to_sdf.py               Convert top-2000 SMILES to 2D SDF
│
├── 04_Ligand_preparation/          3D conformer generation and ROC ligand prep
│   ├── ligand_preparation_pipeline.ipynb    Top-2000 → 3D SDF (ETKDGv3 + MMFF94s)
│   └── roc_ligand_preparation.ipynb         Actives/decoys SMILES → SDF for ROC
│
├── 05_BRD4_DOCK/                   BRD4 docking, ROC validation, RF-Score rescoring
│   ├── 5013_Project_Docking_BRD4.ipynb     Protein prep + ligand prep + cognate redocking
│   ├── run_docking.py              HPC docking script (1955 compounds, Vinardo ex8)
│   ├── submit_docking.pbs          PBS submission for production docking
│   ├── config.txt                  Default Vina configuration
│   ├── submit_all.sh               Auto-submits 6 ROC validation jobs
│   ├── run_roc.pbs                 Single ROC docking job template
│   ├── run_rfscore.pbs             RF-Score-VS rescoring job
│   ├── merge_scores.py             Merges Vina + RF-Score results
│   └── roc_ef_plot.ipynb           Generate ROC/EF curves
│
├── 06_BRD2_DOCK/                   BRD2 docking and ROC validation
│   ├── 5013_Project_Docking_BRD2_v1.ipynb
│   ├── run_docking.py
│   ├── run_roc.pbs
│   └── config.txt
│
├── 07_MD_simulation/               MD simulation (Amber 22)
│   └── (MD scripts and analysis notebooks)
│
├── requirements.txt                Python dependencies
└── README.md                       This file
```

---

## Pipeline Workflow

```
ChEMBL bioactivity data (BRD4: 7,634 compounds; BRD2: 806 compounds)
        ↓
01–02:  QSAR model training and validation (RF + ECFP4)
        ↓
03:     Virtual screening of ~1M ZINC22 lead-like compounds → top 2,000
        ↓
04:     Ligand preparation (3D conformers, MMFF94s minimisation) → 1,955 compounds
        ↓
05–06:  Molecular docking against BRD4 + BRD2 (Vinardo ex8)
        Protocol validation: redocking + ROC/EF + RF-Score-VS rescoring
        ↓
        Candidate selection by BRD4/BRD2 ΔG differential
        ↓
07:     MD simulation (100 ns × 3 replicates) + MM-GBSA
        ↓
        Final candidates: LIG_1629, LIG_1610, LIG_0061
```

---

## Installation

### Python environment

```bash
git clone https://github.com/run209/BRD4-Selective-Inhibitor-VS.git
cd BRD4-Selective-Inhibitor-VS
pip install -r requirements.txt
```

### External tools required

| Tool | Used in | Notes |
|---|---|---|
| AutoDock Vina (≥1.2.3) | 05, 06 | Python bindings via `vina` package |
| MGLTools | 05, 06 | `prepare_receptor4.py` for PDB → PDBQT conversion |
| OpenBabel | 04, 05, 06 | Ligand format conversion |
| RF-Score-VS | 05 | Random Forest rescoring of docking poses |
| AmberTools 23 / Amber 22 | 07 | MD simulation |

### HPC environment (for steps 03, 05, 06)

The virtual screening (~1M compounds) and large-scale docking (1,955 compounds × 2 targets) were run on a PBS-managed HPC cluster. Single-core jobs with 8 GB RAM are sufficient; sample PBS scripts are provided in each folder.

---

## Usage

### 1. Train QSAR models (01, 02)

Open the notebooks in Jupyter or Google Colab:

```bash
jupyter notebook 01_BRD4_QSAR_MODEL/QSAR_BRD4.ipynb
```

These notebooks handle ChEMBL data retrieval, structure curation, scaffold-based splitting, model training (RF/SVR/XGBoost/GB on ECFP4 and 2D descriptors), and full validation (Golbraikh-Tropsha, Y-randomisation, bootstrap CI).

### 2. Virtual screening of ZINC22 (03)

On the HPC cluster:

```bash
cd 03_BRD4_VS_1M/
qsub run_vs_brd4.sh
```

This runs `hpc_mega_screen.py` to apply the trained BRD4 RF model to ~1M compounds, outputting `BRD4_1M_Top2000.csv` (top 2,000 by predicted pIC50, with applicability domain reliability flags).

Then convert top-2000 SMILES to SDF:

```bash
python csv_to_sdf.py BRD4_1M_Top2000.csv BRD4_library.sdf
```

### 3. Ligand preparation (04)

Open `ligand_preparation_pipeline.ipynb` in Colab to generate 3D conformers (ETKDGv3, 10 conformers per molecule, MMFF94s minimisation, lowest-energy conformer retained). Output: `BRD4_library_clean.sdf` (1,955 successfully prepared compounds).

For ROC validation actives/decoys, run `roc_ligand_preparation.ipynb` first to convert the SMILES to SDF.

### 4. Docking (05, 06)

**Protein and ligand preparation, plus cognate redocking validation** — run the Jupyter notebook (`5013_Project_Docking_BRD4.ipynb`) in Google Colab. The notebook handles full receptor preparation (PDBFixer → pdb2pqr/PropKa → OpenMM minimisation → MGLTools PDBQT conversion) and generates `protein.pdbqt`.

**Production docking on HPC:**

```bash
cd 05_BRD4_DOCK/
# Upload BRD4_library_clean.sdf and protein.pdbqt to working directory
qsub submit_docking.pbs
```

This runs `run_docking.py` to dock all 1,955 compounds against BRD4 (Vinardo, exhaustiveness=8, seed=42), outputting `summary.csv` with affinity scores and physicochemical properties.

**ROC/EF validation (six parameter combinations):**

```bash
bash submit_all.sh   # auto-submits 6 jobs: vina/vinardo × ex 8/16/32
```

After completion, RF-Score-VS rescoring:

```bash
qsub run_rfscore.pbs
python merge_scores.py    # combines Vina + RF-Score outputs
```

Then plot ROC/EF curves using `roc_ef_plot.ipynb`.

For BRD2, follow the same workflow in `06_BRD2_DOCK/`. Note that BRD2 used different grid box parameters: centred at (13.72, −2.63, 3.37) Å, dimensions 21.70 × 21.77 × 16.59 Å.

### 5. MD simulation (07)

See `07_MD_simulation/` for Amber 22 input files, equilibration scripts, and trajectory analysis notebooks (RMSD, RMSF, H-bond persistence, MM-GBSA).

---

## Key Parameters

### QSAR (BRD4)

- **Features:** ECFP4 (radius 2, 2048 bits)
- **Algorithm:** Random Forest
- **Optimal hyperparameters:** n_estimators=100, max_depth=None, max_features=sqrt, min_samples_split=2, min_samples_leaf=1
- **Test R² = 0.779** (95% CI: 0.748–0.808); Y-randomisation z = 61.26
- **Applicability domain:** kNN Tanimoto similarity (k=5, threshold=0.3)

### Docking

| Target | PDB | Box centre (Å) | Box size (Å) | Scoring | Exhaustiveness |
|---|---|---|---|---|---|
| BRD4 | 3MXF (1.60 Å) | (28.75, 15.83, −2.34) | 19.29 × 19.84 × 20.49 | Vinardo | 8 |
| BRD2 | 7LAH | (13.72, −2.63, 3.37) | 21.70 × 21.77 × 16.59 | Vinardo | 8 |

- **Cognate redocking RMSD:** BRD4 = 0.32 Å; BRD2 = 1.53 Å
- **ROC/EF validation:** BRD4 AUC = 0.864 (RF-Score-VS rescoring); BRD2 AUC = 0.942 (Vina)

---

## Final Candidates

| Compound | Formula | MW (Da) | LogP | TPSA (Å²) | BRD4 (kcal/mol) | BRD2 (kcal/mol) | ΔG (kcal/mol) |
|---|---|---|---|---|---|---|---|
| LIG_1629 | C₁₇H₁₈FN₃O₂S | 347.1 | 3.06 | 88.8 | −8.237 | −5.828 | −2.409 |
| LIG_1610 | C₁₈H₁₇N₃O₃ | 323.1 | 2.74 | 65.7 | −7.373 | −5.589 | −1.784 |
| LIG_0061 | C₁₆H₁₉F₂N₃O₂S | 355.1 | 3.13 | 64.8 | −7.194 | −5.201 | −1.993 |

All candidates satisfy Lipinski's Rule of Five.

---

## Data Availability

Large data files (raw docking results, ZINC22 library, MD trajectories) are not hosted on GitHub due to size constraints. Code-reproducible data can be regenerated by running the notebooks in order (01 → 07).

The protein PDBQT files are not committed to this repository; they can be regenerated by running Step 1 of the docking notebooks (`5013_Project_Docking_BRD*.ipynb`). The receptor file `protein.pdbqt` used in `05_BRD4_DOCK/` and `09_apo_protein_clean_FINAL.pdbqt` used in ROC validation are the same file; the latter is the named output from the protein preparation pipeline, renamed to `protein.pdbqt` for use with `run_docking.py`.

---

## Notes on Reproducibility

- **Random seeds:** All Vina docking runs use `seed = 42` (per-ligand seed = 42 + global_idx for parallel chunking)
- **Scaffold splitting:** QSAR train/val/test partitions verified to have zero scaffold overlap
- **Y-randomisation:** 30 iterations following Tropsha (2010); seed = `iteration_index`
- **Bootstrap CI:** 1,000 resamples with stratified sampling

Software versions are documented in `requirements.txt` and `Supporting Information S1` of the manuscript.

---

## Citation

If you use this pipeline, please cite the accompanying manuscript:

> [Author names]. Selective BRD4 Inhibitor Discovery via Integrated Virtual Screening. PHM5013 Precision Drug Discovery and Pharmacogenomics, 2026.

---

## Contact

For questions about this pipeline, please contact the corresponding author or open an issue on this repository.
