#!/bin/bash

# ==========================================
# Preprocessing Scripts & One-Linear Command
# ==========================================

# Generate mol2 files with Gasteiger charges using antechamber
antechamber -i BRD4_1629_addH.sdf -fi sdf -o BRD4_1629.mol2 -fo mol2 -c gas -at gaff2


# Generate frcmod files using parmchk2
parmchk2 -i BRD4_1629.mol2 -f mol2 -o BRD4_1629.frcmod


# Generate final topological and coordinate files
tleap -f leap.in

# Run MD simulation
qsub run_amber.pbs

# Trajectory Concatenation, Mergigng, Topology Conversion
cpptraj -i concat_traj.in
cpptraj -i merge.in
cpptraj -i top_conv.in

# Analyze trajectory
cpptraj -i analysis.in

# MMPBSA
qsub run_mmpbsa.pbs

