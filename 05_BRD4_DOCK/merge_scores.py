import os
import io

print("1. Scanning PDBQT files to build name mapping...")
mapping = {}
for f in os.listdir('.'):
    if f.endswith('_out.pdbqt'):
        base_name = f.replace('_out.pdbqt', '')
        internal_name = None
        with io.open(f, 'r', encoding='utf-8', errors='ignore') as file:
            for line in file:
                if 'mol_' in line:
                    for w in line.split():
                        if 'mol_' in w:
                            internal_name = w.strip()
                            break
                if internal_name:
                    break
        if internal_name:
            mapping[base_name] = internal_name

print("2. Reading Vina scores from summary.csv...")
vina_scores = {}
with open('summary.csv', 'r') as f:
    lines = f.readlines()
    for line in lines[1:]:
        parts = line.strip().split(',')
        if len(parts) >= 3:
            lname = parts[0]
            # Skip entries with '_out' suffix (artefacts)
            if not lname.endswith('_out'):
                vina_scores[lname] = {'Is_active': parts[1], 'Affinity': parts[2]}

print("3. Reading RF-Score-VS scores from rfscore_summary.csv...")
rf_scores = {}
with open('rfscore_summary.csv', 'r') as f:
    lines = f.readlines()
    for line in lines[1:]:
        parts = line.strip().split(',')
        if len(parts) >= 3:
            mol_name = parts[1]
            try:
                score = float(parts[2])
                # Keep the highest score across poses for the same molecule
                if mol_name not in rf_scores or score > rf_scores[mol_name]:
                    rf_scores[mol_name] = score
            except:
                pass

print("4. Merging into final_master_scores.csv...")
with open('final_master_scores.csv', 'w') as f:
    f.write("Ligand_name,name,Is_active,Affinity,RF_Score\n")
    for lname, v_data in vina_scores.items():
        if lname in mapping:
            mname = mapping[lname]
            if mname in rf_scores:
                rf_val = rf_scores[mname]
                f.write("{},{},{},{},{}\n".format(
                    lname, mname, v_data['Is_active'], v_data['Affinity'], rf_val
                ))

print("Done. final_master_scores.csv generated successfully.")
