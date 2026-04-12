# =============================================================================
#  AutoDock Vina -- Docking Pipeline for HPC
# =============================================================================

import argparse, os, re, csv, time, zipfile
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors
from openbabel import openbabel as ob
from vina import Vina

RDLogger.DisableLog("rdApp.*")

# ---------------------------------------------------------------------------
#  Fixed parameters (from cognate docking validation)
# ---------------------------------------------------------------------------

CENTER      = [13.72, -2.63, 3.37]
SIZE        = [21.70, 21.77, 16.59]
SCORING_FN  = "vinardo"
EXHAUSTIVE  = 8
N_POSES     = 1
SEED        = 42
SKIP_FAILED = True

# ---------------------------------------------------------------------------
#  Argument parser
# ---------------------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("--protein",  required=True, help="protein.pdbqt")
parser.add_argument("--ligands",  required=True, help="ligands.sdf")
parser.add_argument("--outdir",   default="docking_results")
parser.add_argument("--start",    type=int, default=0,   help="start index")
parser.add_argument("--end",      type=int, default=None, help="end index")
args = parser.parse_args()

OUTDIR = Path(args.outdir)
OUTDIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
#  Helper functions
# ---------------------------------------------------------------------------

def sdf_to_pdbqt(sdf_path, pdbqt_path):
    obc = ob.OBConversion()
    obc.SetInAndOutFormats("sdf", "pdbqt")
    mol = ob.OBMol()
    if not obc.ReadFile(mol, str(sdf_path)):
        raise RuntimeError(f"Cannot read {sdf_path}")
    cm = ob.OBChargeModel.FindType("gasteiger")
    if cm:
        cm.ComputeCharges(mol)
    obc.WriteFile(mol, str(pdbqt_path))

def pdbqt_to_sdf(pdbqt_path, sdf_path):
    obc = ob.OBConversion()
    obc.SetInAndOutFormats("pdbqt", "sdf")
    mol = ob.OBMol()
    obc.ReadFile(mol, str(pdbqt_path))
    obc.WriteFile(mol, str(sdf_path))

# ---------------------------------------------------------------------------
#  Load ligands
# ---------------------------------------------------------------------------

print(f"Loading ligands from {args.ligands} ...")
mol_list = []
for idx, mol in enumerate(Chem.SDMolSupplier(args.ligands, removeHs=False, sanitize=True)):
    if mol is None:
        continue
    raw  = mol.GetPropsAsDict().get("_Name", "").strip()
    name = re.sub(r"[^\w\-]", "_", raw) if raw else f"LIG_{idx+1:04d}"
    mol_list.append((name, mol))

# Apply index range
end = args.end if args.end else len(mol_list)
mol_list = mol_list[args.start:end]
print(f"Docking {len(mol_list)} ligands (index {args.start} to {end})")

# ---------------------------------------------------------------------------
#  Docking loop
# ---------------------------------------------------------------------------

summary_rows = []
t0_total = time.time()

for lig_idx, (lig_name, mol_raw) in enumerate(mol_list):
    global_idx = args.start + lig_idx
    print(f"\n[{global_idx+1}] {lig_name}")

    try:
        mw      = Descriptors.ExactMolWt(mol_raw)
        logp    = Descriptors.MolLogP(mol_raw)
        hbd     = rdMolDescriptors.CalcNumHBD(mol_raw)
        hba     = rdMolDescriptors.CalcNumHBA(mol_raw)
        rotb    = rdMolDescriptors.CalcNumRotatableBonds(mol_raw)
        tpsa    = Descriptors.TPSA(mol_raw)
        formula = rdMolDescriptors.CalcMolFormula(mol_raw)
        ro5_ok  = sum([mw<=500, logp<=5, hbd<=5, hba<=10]) >= 3

        lig_dir   = OUTDIR / f"ligand_{global_idx+1:04d}_{lig_name}"
        lig_dir.mkdir(exist_ok=True)
        lig_sdf   = lig_dir / "ligand.sdf"
        lig_pdbqt = lig_dir / "ligand.pdbqt"
        out_pdbqt = lig_dir / "docked.pdbqt"
        out_sdf   = lig_dir / "docked.sdf"

        w = Chem.SDWriter(str(lig_sdf))
        w.write(mol_raw)
        w.close()

        sdf_to_pdbqt(lig_sdf, lig_pdbqt)

        v = Vina(sf_name=SCORING_FN, seed=SEED+global_idx, verbosity=0)
        v.set_receptor(args.protein)
        v.set_ligand_from_file(str(lig_pdbqt))
        v.compute_vina_maps(center=CENTER, box_size=SIZE)

        t0 = time.time()
        v.dock(exhaustiveness=EXHAUSTIVE, n_poses=N_POSES)
        elapsed = time.time() - t0

        energies = v.energies(n_poses=N_POSES)
        best_aff = energies[0][0]
        print(f"  Affinity = {best_aff:.3f} kcal/mol  ({elapsed:.1f}s)")

        v.write_poses(str(out_pdbqt), n_poses=N_POSES, overwrite=True)
        pdbqt_to_sdf(out_pdbqt, out_sdf)

        summary_rows.append({
            "Ligand_name": lig_name,
            "Formula": formula, "MW_Da": round(mw,3), "LogP": round(logp,3),
            "HBD": hbd, "HBA": hba, "RotB": rotb, "TPSA_A2": round(tpsa,2),
            "Lipinski_Ro5": "PASS" if ro5_ok else "FAIL",
            "Best_affinity_kcal_mol": round(best_aff,4),
            "Runtime_s": round(elapsed,2), "Status": "SUCCESS",
        })

    except Exception as err:
        print(f"  FAILED: {err}")
        summary_rows.append({
            "Ligand_name": lig_name, "Formula": "--", "MW_Da": "--",
            "LogP": "--", "HBD": "--", "HBA": "--", "RotB": "--",
            "TPSA_A2": "--", "Lipinski_Ro5": "--",
            "Best_affinity_kcal_mol": "--", "Runtime_s": "--",
            "Status": f"FAILED: {err}",
        })
        if not SKIP_FAILED:
            raise

# ---------------------------------------------------------------------------
#  Save results
# ---------------------------------------------------------------------------

successful = [r for r in summary_rows if r["Status"] == "SUCCESS"]
successful.sort(key=lambda r: r["Best_affinity_kcal_mol"])
for rank, row in enumerate(successful, 1):
    row["Rank"] = rank

csv_path = OUTDIR / f"summary_{args.start}_{end}.csv"
fields = ["Rank", "Ligand_name", "Formula", "MW_Da", "LogP", "HBD", "HBA",
          "RotB", "TPSA_A2", "Lipinski_Ro5", "Best_affinity_kcal_mol",
          "Runtime_s", "Status"]
with open(csv_path, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=fields)
    writer.writeheader()
    writer.writerows(summary_rows)

total = time.time() - t0_total
print(f"\nDone! {len(successful)}/{len(summary_rows)} success | {total/60:.1f} min")
print(f"Results: {csv_path}")