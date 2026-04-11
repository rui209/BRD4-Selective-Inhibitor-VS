#!/bin/bash
# submit_all.sh
# Submits 6 ROC docking jobs: vina/vinardo x exhaustiveness 8/16/32

WORKDIR=/home/svu/e1569136
SRCDIR=${WORKDIR}/ROC_BRD4_V2
RECEPTOR=${SRCDIR}/09_apo_protein_clean_FINAL.pdbqt

for scoring in vina vinardo; do
    for exh in 8 16 32; do

        JOBNAME="ROC_${scoring}_ex${exh}"
        JOBDIR="${WORKDIR}/${JOBNAME}"

        mkdir -p ${JOBDIR}

        # Copy ligand files and Vina binary
        cp ${SRCDIR}/active_*.pdbqt ${JOBDIR}/
        cp ${SRCDIR}/decoy_*.pdbqt ${JOBDIR}/
        cp ${SRCDIR}/vina ${JOBDIR}/

        # Generate config.txt
        cat > ${JOBDIR}/config.txt << EOF
receptor = ${RECEPTOR}
scoring = ${scoring}
exhaustiveness = ${exh}
center_x = 28.75
center_y = 15.83
center_z = -2.34
size_x = 19.29
size_y = 19.84
size_z = 20.49
EOF

        # Generate PBS script
        cat > ${JOBDIR}/run.pbs << EOF
#!/bin/bash
#PBS -N ${JOBNAME}
#PBS -l nodes=1:ppn=16
#PBS -l walltime=48:00:00
#PBS -l mem=32gb
#PBS -j oe

cd ${JOBDIR}

echo "Ligand_name,Is_active,Affinity" > summary.csv

for ligand in active_*.pdbqt decoy_*.pdbqt; do
    [ -e "\$ligand" ] || continue
    name=\$(basename "\$ligand" .pdbqt)
    if [[ "\$name" == active_* ]]; then
        is_active=1
    else
        is_active=0
    fi
    ${JOBDIR}/vina \\
        --config ${JOBDIR}/config.txt \\
        --ligand "\$ligand" \\
        --out "\${name}_out.pdbqt" > "vina_\${name}.log" 2>&1
    score=\$(awk '/^ *1 +[-0-9]/{print \$2; exit}' "vina_\${name}.log")
    if [ -z "\$score" ]; then score=0.0; fi
    echo "\$name,\$is_active,\$score" >> summary.csv
done

echo "Done: ${JOBNAME}"
EOF

        sed -i 's/\r//' ${JOBDIR}/run.pbs
        qsub ${JOBDIR}/run.pbs
        echo "Submitted: ${JOBNAME}"

    done
done
