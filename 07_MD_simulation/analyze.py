from google.colab import drive
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob

drive.mount('/content/drive')

base_path = '/content/drive/MyDrive/5013_project_BRD4/07_MD_simulation/final'

if not os.path.exists(base_path):
    print(f"warning：{base_path} do not exist")
else:
    print("successfully connected")

def plot_rmsd_group(system_prefix):
    plt.figure(figsize=(10, 6))
    all_reps = []

    folders = sorted(glob.glob(os.path.join(base_path, f"{system_prefix}_*")))
    if not folders:
        print(f"Error： {base_path} not found {system_prefix}_*")
        return

    processed_reps = []
    for i, folder in enumerate(folders):
        file_path = os.path.join(folder, 'rmsd_protein.dat')
        if os.path.exists(file_path):
            try:
                df = pd.read_csv(file_path, sep=r'\s+', header=None, comment='#', on_bad_lines='skip')

                if df.shape[1] >= 2:
                    data = pd.to_numeric(df.iloc[:, 1], errors='coerce').dropna().values
                    processed_reps.append(data)
                    plt.plot(range(len(data)), data, alpha=0.3, label=f"Rep {i+1}")
                else:
                    print(f"warning：{file_path} not enough")
            except Exception as e:
                print(f"cannot read {file_path}: {e}")

    if processed_reps:
        min_len = min(len(r) for r in processed_reps)

        truncated_data = np.array([r[:min_len] for r in processed_reps])

        mean_val = np.mean(truncated_data, axis=0)
        std_val = np.std(truncated_data, axis=0)
        frames = np.arange(min_len)

        color_reps = 'skyblue'
        color_mean = 'tab:blue'
        color_std = 'gainsboro'

        plt.fill_between(frames, mean_val - std_val, mean_val + std_val,
                         color=color_std, alpha=0.5, zorder=1, label='Std. Dev.')
        for i, data in enumerate(processed_reps[:3]):
            plt.plot(frames, data[:min_len], color=color_reps, alpha=0.3, linewidth=0.8, zorder=2)

        plt.plot(frames, mean_val, color=color_mean, linewidth=2.5, zorder=3, label='Average')

        plt.title(f'RMSD Stability & Error Analysis: {system_prefix}', fontsize=16, fontweight='bold')
        plt.xlabel('Time (Frames)', fontsize=14)
        plt.ylabel('RMSD (Å)', fontsize=14)
        plt.legend(loc='upper right', fontsize=11, frameon=False, title_fontsize=12)
        plt.grid(axis='both', linestyle='--', alpha=0.5)

        sns.despine()

        out_name = f'{system_prefix}_rmsd_pro_color.png'
        save_path = os.path.join(base_path, out_name)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"figure saved: {out_name}")
        plt.show()

plot_rmsd_group('BRD4_1610')

systems = ['BRD4_1629', 'BRD4_1610', 'BRD4_0061', 'BRD2_1629', 'BRD2_1610', 'BRD2_0061']
summary_data = []

for sys in systems:
    rep_vals = []
    folders = sorted(glob.glob(os.path.join(base_path, f"{sys}_*")))

    for folder in folders:
        mmpbsa_file = os.path.join(folder, 'mmpbsa_results.dat')
        if os.path.exists(mmpbsa_file):
            with open(mmpbsa_file, 'r') as f:
                for line in f:
                    if "DELTA TOTAL" in line:
                        val = float(line.split()[2])
                        rep_vals.append(val)
                        break

    if rep_vals:
        mean_g = np.mean(rep_vals)
        std_g = np.std(rep_vals)
        summary_data.append({
            'System': sys,
            'Mean_dG': mean_g,
            'Std_Dev': std_g,
            'Reps': len(rep_vals)
        })

df_final = pd.DataFrame(summary_data)
print(df_final)

df_final.to_csv(os.path.join(base_path, 'mmpbsa_summary_table.csv'), index=False)

def plot_mmpbsa_clean(df):
    df_plot = df.copy()
    df_plot['Protein'] = df_plot['System'].str[:4]
    df_plot['Conformer'] = df_plot['System'].str[5:]

    sns.set_style("white")
    plt.figure(figsize=(10, 6))

    ax = sns.barplot(
        data=df_plot,
        x='Conformer',
        y='Mean_dG',
        hue='Protein',
        palette="muted",
        edgecolor='black',
        linewidth=1.2,
        capsize=.1,
        errwidth=1.5
    )


    x_coords = [p.get_x() + p.get_width() / 2. for p in ax.patches if p.get_height() != 0]

    df_resorted = pd.concat([
        df_plot[df_plot['Protein'] == 'BRD4'].sort_values('Conformer'),
        df_plot[df_plot['Protein'] == 'BRD2'].sort_values('Conformer')
    ])

    plt.errorbar(
        x=x_coords,
        y=df_resorted['Mean_dG'],
        yerr=df_resorted['Std_Dev'],
        fmt='none',
        c='black',
        capsize=4,
        elinewidth=1.2
    )

    plt.axhline(0, color='black', linewidth=1)
    plt.ylabel(r'$\Delta G_{bind}$ (kcal/mol)', fontsize=14)
    plt.xlabel('Ligand Conformer', fontsize=14)
    plt.title('Binding Free Energy with SD (Replicates n=3)', fontsize=15, fontweight='bold')

    plt.ylim(df_plot['Mean_dG'].min() - 8, 2)

    plt.legend(title='Target Protein', frameon=False, loc='lower right')
    sns.despine()
    plt.grid(axis='y', linestyle=':', alpha=0.5)

    save_path = os.path.join(base_path, 'mmpbsa_final_report.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"✅ saved {save_path}")
    plt.show()

plot_mmpbsa_clean(df_final)

def plot_hb_count_heatmap(systems):
    all_hb_counts = {}

    for sys in systems:
        rep_vals = []
        folders = sorted(glob.glob(os.path.join(base_path, f"{sys}_*")))

        for folder in folders:
            file_path = os.path.join(folder, 'nhb_prot2lig.dat')
            if os.path.exists(file_path):
                df = pd.read_csv(file_path, sep=r'\s+')
                avg_count = df.iloc[:, 1].mean()
                rep_vals.append(avg_count)

        if rep_vals:
            all_hb_counts[sys] = np.mean(rep_vals)

    res_df = pd.DataFrame([all_hb_counts])

    plt.figure(figsize=(10, 3))
    sns.heatmap(res_df, annot=True, fmt=".2f", cmap="YlGnBu", cbar_kws={'label': 'Avg. HB Count'})
    plt.title('Average Hydrogen Bond Count (n=3)')
    plt.xlabel('System')
    plt.ylabel('')
    save_path = os.path.join(base_path, 'H-bonds comparison.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()

all_systems = ['BRD4_1629', 'BRD4_1610', 'BRD4_0061', 'BRD2_1629', 'BRD2_1610', 'BRD2_0061']
plot_hb_count_heatmap(all_systems)

def plot_rmsf_absolute_clean(sys1, sys2):
    sns.set_style("white")
    plt.figure(figsize=(12, 6))

    def get_clean_data(name):
        folders = sorted(glob.glob(os.path.join(base_path, f"{name}_*")))
        processed_reps = []
        min_len = float('inf')
        res_indices = None

        for f in folders:
            path = os.path.join(f, 'rmsf.dat')
            if os.path.exists(path):
                df = pd.read_csv(path, sep=r'\s+', header=None, comment='#')


                res_col = pd.to_numeric(df.iloc[:, 0], errors='coerce')
                val_col = pd.to_numeric(df.iloc[:, 1], errors='coerce')

                clean_df = pd.concat([res_col, val_col], axis=1).dropna()

                val = clean_df.iloc[:, 1].values
                processed_reps.append(val)
                min_len = min(min_len, len(val))
                res_indices = clean_df.iloc[:, 0].values

        if not processed_reps:
            print(f"⚠️ not found: {name}")
            return None, None, None

        aligned_reps = [r[:min_len] for r in processed_reps]
        aligned_res = res_indices[:min_len]

        mean_val = np.mean(aligned_reps, axis=0)
        std_val = np.std(aligned_reps, axis=0)
        return aligned_res, mean_val, std_val

    res1, mean1, std1 = get_clean_data(sys1)
    res2, mean2, std2 = get_clean_data(sys2)

    if res1 is not None:
        plt.plot(res1, mean1, label=sys1, color='tab:blue', lw=2)
        plt.fill_between(res1, mean1-std1, mean1+std1, color='tab:blue', alpha=0.2)
    if res2 is not None:
        plt.plot(res2, mean2, label=sys2, color='tab:orange', lw=2)
        plt.fill_between(res2, mean2-std2, mean2+std2, color='tab:orange', alpha=0.2)

    plt.title('Protein Residue Fluctuation (RMSF)', fontsize=15, fontweight='bold')
    plt.xlabel('Residue Index', fontsize=12)
    plt.ylabel('RMSF (Å)', fontsize=12)
    plt.legend(frameon=False)
    sns.despine()
    save_path = os.path.join(base_path, 'RMSF(BRD4_1629vsBRD2_1629).png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()

plot_rmsf_absolute_clean('BRD4_1629', 'BRD2_1629')

def plot_rmsf_single_protein_comparison(protein_prefix, ligands):
    sns.set_style("white")
    plt.figure(figsize=(12, 6))

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

    for idx, lig in enumerate(ligands):
        sys_name = f"{protein_prefix}_{lig}"
        folders = sorted(glob.glob(os.path.join(base_path, f"{sys_name}_*")))

        processed_reps = []
        res_indices = None

        for f in folders:
            path = os.path.join(f, 'rmsf.dat')
            if os.path.exists(path):
                df = pd.read_csv(path, sep=r'\s+', header=None, comment='#')
                res_col = pd.to_numeric(df.iloc[:, 0], errors='coerce')
                val_col = pd.to_numeric(df.iloc[:, 1], errors='coerce')
                clean_df = pd.concat([res_col, val_col], axis=1).dropna()

                processed_reps.append(clean_df.iloc[:, 1].values)
                res_indices = clean_df.iloc[:, 0].values

        if processed_reps:
            min_len = min(len(r) for r in processed_reps)
            aligned_reps = [r[:min_len] for r in processed_reps]
            mean_val = np.mean(aligned_reps, axis=0)
            std_val = np.std(aligned_reps, axis=0)

            label_name = f"LIG {lig}"
            plt.plot(res_indices[:min_len], mean_val, label=label_name, color=colors[idx], lw=2)
            plt.fill_between(res_indices[:min_len], mean_val - std_val, mean_val + std_val,
                             color=colors[idx], alpha=0.15)
        else:
            print(f"⚠️ not found {sys_name}")

    plt.title(f'Protein Residue Fluctuation (RMSF) on {protein_prefix}', fontsize=16, fontweight='bold')
    plt.xlabel('Residue Index', fontsize=13)
    plt.ylabel('RMSF (Å)', fontsize=13)

    plt.legend(frameon=False, fontsize=11)
    plt.grid(axis='y', linestyle=':', alpha=0.5)
    sns.despine()

    save_path = os.path.join(base_path, f'rmsf_{protein_prefix}_ligand_comparison.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()


plot_rmsf_single_protein_comparison('BRD4', ['1629', '1610', '0061'])

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import glob

def plot_final_rmsd_matrix(proteins, ligands, base_path):
    sns.set_style("white")
    fig, axes = plt.subplots(len(proteins), len(ligands), figsize=(18, 10), sharey=True)

    color_p = '#1f77b4'
    color_l = '#ff7f0e'

    for p_idx, protein in enumerate(proteins):
        for l_idx, ligand in enumerate(ligands):
            ax = axes[p_idx, l_idx]
            sys_name = f"{protein}_{ligand}"
            folders = sorted(glob.glob(os.path.join(base_path, f"{sys_name}_*")))

            p_reps = []
            l_reps = []
            min_len = float('inf')

            for f in folders:
                p_path = os.path.join(f, 'rmsd_protein.dat')
                l_path = os.path.join(f, 'rmsd_ligand.dat')

                if os.path.exists(p_path) and os.path.exists(l_path):
                    try:
                        p_df = pd.read_csv(p_path, sep=r'\s+', header=None, comment='#')
                        l_df = pd.read_csv(l_path, sep=r'\s+', header=None, comment='#')

                        p_vals = pd.to_numeric(p_df.iloc[:, 1], errors='coerce').dropna().values
                        l_vals = pd.to_numeric(l_df.iloc[:, 1], errors='coerce').dropna().values

                        curr_len = min(len(p_vals), len(l_vals))
                        if curr_len > 0:
                            p_reps.append(p_vals[:curr_len])
                            l_reps.append(l_vals[:curr_len])
                            min_len = min(min_len, curr_len)
                    except:
                        continue

            if p_reps and l_reps:
                time_axis = np.arange(min_len) * 0.1

                p_mean = np.mean([r[:min_len] for r in p_reps], axis=0)
                p_std = np.std([r[:min_len] for r in p_reps], axis=0)

                l_mean = np.mean([r[:min_len] for r in l_reps], axis=0)
                l_std = np.std([r[:min_len] for r in l_reps], axis=0)

                ax.plot(time_axis, p_mean, color=color_p, label=f'Protein (mean {p_mean.mean():.2f} Å)', lw=1.5)
                ax.fill_between(time_axis, p_mean - p_std, p_mean + p_std, color=color_p, alpha=0.15)

                ax.plot(time_axis, l_mean, color=color_l, label=f'Ligand (mean {l_mean.mean():.2f} Å)', lw=1.5)
                ax.fill_between(time_axis, l_mean - l_std, l_mean + l_std, color=color_l, alpha=0.15)

            ax.set_title(f"{sys_name}", fontsize=14, fontweight='bold')
            ax.legend(fontsize='small', loc='upper left', frameon=True)
            if p_idx == 1: ax.set_xlabel('Time (ns)', fontsize=12)
            if l_idx == 0: ax.set_ylabel('RMSD (Å)', fontsize=12)
            ax.set_ylim(0, 6)

    plt.suptitle('Structural Stability (RMSD) Comparison: Protein vs Ligand', fontsize=18, fontweight='bold', y=1.02)
    plt.tight_layout()

    save_path = os.path.join(base_path, 'rmsd_final_matrix_combined.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"✅ saved {save_path}")
    plt.show()

plot_final_rmsd_matrix(['BRD4', 'BRD2'], ['1629', '1610', '0061'], base_path)

import pandas as pd
import os
import glob

def extract_mmpbsa_components(proteins, ligands, base_path):
    results = []

    for protein in proteins:
        for ligand in ligands:
            sys_name = f"{protein}_{ligand}"
            folders = glob.glob(os.path.join(base_path, f"{sys_name}_*"))

            rep_data = {'VDW': [], 'ELE': [], 'SOLV': []}

            for f in folders:
                file_path = os.path.join(f, 'mmpbsa_results.dat')

                if os.path.exists(file_path):
                    with open(file_path, 'r') as file:
                        lines = file.readlines()
                        for line in lines:
                            if 'VDW' in line and 'Δ' not in line:
                                try:
                                    val = float(line.split()[1])
                                    rep_data['VDW'].append(val)
                                except: continue
                            if 'EEL' in line:
                                try:
                                    val = float(line.split()[1])
                                    rep_data['ELE'].append(val)
                                except: continue
                            if 'EGB' in line or 'EPB' in line:
                                try:
                                    egb = float(line.split()[1])
                                    rep_data['SOLV'].append(egb)
                                except: continue

            if rep_data['VDW']:
                for comp in ['VDW', 'ELE', 'SOLV']:
                    results.append({
                        'System': sys_name,
                        'Component': comp,
                        'Energy': sum(rep_data[comp]) / len(rep_data[comp])
                    })

    return pd.DataFrame(results)

proteins = ['BRD4', 'BRD2']
ligands = ['1629', '1610', '0061']

df_energy = extract_mmpbsa_components(proteins, ligands, base_path)

print(df_energy.head(10))

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob

def get_cleaned_energy_data(proteins, ligands, base_path):
    final_data = []

    for protein in proteins:
        for ligand in ligands:
            sys_name = f"{protein}_{ligand}"
            search_pattern = os.path.join(base_path, f"{sys_name}*")
            folders = [f for f in glob.glob(search_pattern) if os.path.isdir(f)]

            print(f"Reading {sys_name}: {len(folders)} folders found.")

            vdw_list, net_ele_list, surf_list = [], [], []

            for f in folders:
                file_path = os.path.join(f, 'mmpbsa_results.dat')
                if not os.path.exists(file_path): continue

                raw = {}
                is_delta_section = False
                with open(file_path, 'r') as file:
                    for line in file:
                        if "Differences (Complex - Receptor - Ligand):" in line:
                            is_delta_section = True
                            continue

                        if is_delta_section:
                            parts = line.split()
                            if len(parts) < 2: continue

                            key = parts[0]
                            try:
                                val = float(parts[1])
                                if key == 'VDWAALS': raw['VDW'] = val
                                if key == 'EEL': raw['EEL'] = val
                                if key == 'EGB': raw['EGB'] = val
                                if key == 'ESURF': raw['ESURF'] = val
                            except: continue

                if all(k in raw for k in ['VDW', 'EEL', 'EGB']):
                    vdw_list.append(raw['VDW'])
                    net_ele_list.append(raw['EEL'] + raw['EGB'])
                    surf_list.append(raw.get('ESURF', 0))

            if vdw_list:
                final_data.append({'System': sys_name, 'Component': 'Van der Waals', 'Energy': np.mean(vdw_list)})
                final_data.append({'System': sys_name, 'Component': 'Net Electrostatic', 'Energy': np.mean(net_ele_list)})
                final_data.append({'System': sys_name, 'Component': 'Non-polar Solvation', 'Energy': np.mean(surf_list)})

    return pd.DataFrame(final_data)

def plot_energy_decomposition(df, save_path):
    if df.empty:
        print("🚨 failed")
        return

    pivot_df = df.pivot(index='System', columns='Component', values='Energy')
    cols = ['Van der Waals', 'Net Electrostatic', 'Non-polar Solvation']
    pivot_df = pivot_df.reindex(columns=cols)

    plt.figure(figsize=(12, 7))
    sns.set_style("white")

    systems = pivot_df.index
    vdw = pivot_df['Van der Waals'].values
    ele = pivot_df['Net Electrostatic'].values
    surf = pivot_df['Non-polar Solvation'].values
    total = vdw + ele + surf

    plt.bar(systems, vdw, label='Van der Waals', color='#3498db', alpha=0.8)
    plt.bar(systems, ele, bottom=vdw, label='Net Electrostatic (EEL+EGB)', color='#2ecc71', alpha=0.8)
    plt.bar(systems, surf, bottom=vdw+ele, label='Non-polar Solvation', color='#95a5a6', alpha=0.8)

    plt.scatter(systems, total, color='black', marker='D', s=100, label=r'Total $\Delta G_{bind}$', zorder=5)

    for i, val in enumerate(total):
        plt.text(i, val - 1.5, f'{val:.2f}', ha='center', fontweight='bold', fontsize=10)

    plt.axhline(0, color='black', lw=1)
    plt.ylabel('Energy (kcal/mol)', fontsize=12)
    plt.title('MMPBSA Energy Decomposition (Corrected Delta G)', fontsize=14, fontweight='bold')
    plt.xticks(rotation=30, ha='right')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, 'mmpbsa_fixed_plot.png'), dpi=300)
    plt.show()


df_energy = get_cleaned_energy_data(['BRD4', 'BRD2'], ['1629', '1610', '0061'], base_path)
plot_energy_decomposition(df_energy, base_path)

def get_persistence_data_perfect(proteins, ligands, base_path):
    all_data = []
    for protein in proteins:
        for ligand in ligands:
            sys_name = f"{protein}_{ligand}"
            search_pattern = os.path.join(base_path, f"{sys_name}_*")
            potential_folders = sorted(glob.glob(search_pattern))

            folders = [f for f in potential_folders if os.path.isdir(f) and f[-1].isdigit()]

            print(f"Reading {sys_name}: Found {len(folders)} valid replicate folders.")

            res_dict = {}
            for f in folders:
                for suffix in ['prot2lig', 'lig2prot']:
                    file_path = os.path.join(f, f'avghb_{suffix}.dat')
                    if os.path.exists(file_path):
                        with open(file_path, 'r') as file:
                            for line in file:
                                if line.startswith('#') or not line.strip(): continue
                                parts = line.split()
                                if len(parts) < 5: continue

                                acc, don = parts[0], parts[2]
                                raw_res = acc.split('@')[0] if 'MOL' not in acc else don.split('@')[0]
                                res = re.sub(r'[^a-zA-Z0-9]', '', str(raw_res)).upper() # 清洗残基名

                                try:
                                    frac = float(parts[4]) * 100
                                    if res not in res_dict: res_dict[res] = []
                                    res_dict[res].append(frac)
                                except: continue

            for res, values in res_dict.items():
                if values:
                    avg_p = sum(values) / len(folders)
                    if avg_p > 1.0:
                        all_data.append({'System': sys_name, 'Residue': res, 'Persistence': avg_p})

    return pd.DataFrame(all_data)


df_p_perfect = get_persistence_data_perfect(['BRD4', 'BRD2'], ['1629', '1610', '0061'], base_path)


if not df_p_perfect.empty:
    print("\n check system", df_p_perfect['System'].unique())

    plot_persistence_heatmap_fixed(df_p_perfect, image_save_path)
    plot_persistence_bar_final(df_p_perfect, image_save_path)

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

def plot_unified_residue_charts(df, save_path):
    df_unified = df.copy()

    def apply_offset(row):
        res_name = "".join(filter(str.isalpha, row['Residue']))
        try:
            res_num = int("".join(filter(str.isdigit, row['Residue'])))
        except:
            return row['Residue']


        if 'BRD4' in row['System']:
            unified_num = res_num + 41
        elif 'BRD2' in row['System']:
            unified_num = res_num + 56
        else:
            unified_num = res_num

        return f"{res_name}{unified_num}"

    df_unified['Residue'] = df_unified.apply(apply_offset, axis=1)


    top_res = df_unified.groupby('Residue')['Persistence'].max().nlargest(15).index
    df_plot = df_unified[df_unified['Residue'].isin(top_res)]


    plt.figure(figsize=(12, 8))
    sns.set_style("ticks")
    ax = sns.barplot(
        data=df_plot,
        x='Persistence',
        y='Residue',
        hue='System',
        palette='turbo',
        edgecolor='black'
    )

    plt.title('H-Bond Occupancy (Unified PDB Numbering: BRD4+41, BRD2+56)', fontsize=15, fontweight='bold')
    plt.xlabel('Persistence (%)', fontsize=12)
    plt.ylabel('Residue (Standard PDB Number)', fontsize=12)
    plt.xlim(0, 110)


    for p in ax.patches:
        width = p.get_width()
        if width > 1:
            ax.annotate(f'{width:.1f}%', (width + 1, p.get_y() + p.get_height()/2),
                        va='center', fontsize=9)

    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.savefig(os.path.join(save_path, 'hbond_bar_unified.png'), dpi=300, bbox_inches='tight')
    plt.show()


    plt.figure(figsize=(10, 8))
    pivot_df = df_unified.pivot_table(index='Residue', columns='System', values='Persistence').fillna(0)

    pivot_df = pivot_df.loc[pivot_df.max(axis=1) > 5]

    sns.heatmap(pivot_df, annot=True, fmt=".1f", cmap='YlGnBu', linewidths=.5, cbar_kws={'label': 'Occupancy (%)'})
    plt.title('Hydrogen Bond Heatmap (Standard Numbering)', fontsize=15, fontweight='bold')
    plt.savefig(os.path.join(save_path, 'hbond_heatmap_unified.png'), dpi=300, bbox_inches='tight')
    plt.show()


plot_unified_residue_charts(df_p_perfect, image_save_path)