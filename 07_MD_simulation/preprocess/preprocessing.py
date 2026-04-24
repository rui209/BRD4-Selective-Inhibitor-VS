# convert .pdbqt to .pdb
# run in bash
# obabel -ipdbqt LIG_1629_out.pdbqt -opdb -O LIG_1629_docked.pdb -f 1 -l 1


# add H before proceeding to MD simulation

!pip install rdkit -q

from rdkit import Chem
from rdkit.Chem import AllChem
from google.colab import drive
import os

drive.mount("/content/drive")

INPUT_DIR  = "/content/drive/MyDrive/5013_project_BRD4/docking_results/BRD2/"
OUTPUT_DIR = "/content/drive/MyDrive/5013_project_BRD4/docking_results/BRD2/addH/"
# ────────────────────────────────────────────────────────

os.makedirs(OUTPUT_DIR, exist_ok=True)

PDB_FILES = sorted([
    os.path.join(INPUT_DIR, f)
    for f in os.listdir(INPUT_DIR)
    if f.endswith(".pdb")
])

for i, path in enumerate(PDB_FILES, 1):
    fname = os.path.basename(path)
    stem  = os.path.splitext(fname)[0]
    print(f"[{i}/{len(PDB_FILES)}] {fname} ...", end=" ")

    mol = Chem.MolFromPDBFile(path, removeHs=True, sanitize=True)
    if mol is None:
        print("✗ failed")
        continue

    mol_H = Chem.AddHs(mol, addCoords=True)


    pdb_block = Chem.MolToPDBBlock(mol_H)
    out_path  = os.path.join(OUTPUT_DIR, f"{stem}_addH.pdb")
    with open(out_path, "w") as f:
        f.write(pdb_block)

    print(f"✓ saved {out_path}")

print("\n all done!")



# convert to .sdf

from rdkit import Chem
import os

OUTPUT_DIR = "/content/drive/MyDrive/5013_project_BRD4/docking_results/BRD2/addH/"

for fname in sorted(os.listdir(OUTPUT_DIR)):
    if not fname.endswith("_addH.pdb"): continue
    stem = os.path.splitext(fname)[0]

    mol = Chem.MolFromPDBFile(os.path.join(OUTPUT_DIR, fname), removeHs=False, sanitize=True)
    if mol is None:
        print(f"✗ {fname} failed")
        continue

    sdf_block = Chem.MolToMolBlock(mol)
    with open(os.path.join(OUTPUT_DIR, f"{stem}.sdf"), "w") as f:
        f.write(sdf_block)

    print(f"✓ {stem}.sdf saved")

