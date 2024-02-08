#!/bin/bash

. ${HOME}/miniconda3/etc/profile.d/conda.sh
conda activate py310-strainfilter

coord=.
#coord=MCULE-1005714968.pdbqt
emaxtotalstrain=6.0
emaxtorsion=1.8
# Edit scriptpath to your StrainFIlter installed dirctory
scriptpath=${HOME}/data/scripts/ChemInfTools/apps/strainfilter

pystrainfilter -c ${coord} --emax-total-strain ${emaxtotalstrain} \
    --emax-torsion ${emaxtorsion} --scriptpath ${scriptpath} -s
