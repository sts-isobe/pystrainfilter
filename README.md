# pystrainfilter
python interface of StrainFilter

## Licence
MIT

## Install
```
$ conda create –n py310-strainfilter –c conda-forge
$ conda activate py310-strainfilter
$ conda install –c conda-forge rdkit openbabel

$ git clone https://github.com/docking-org/ChemInfTools.git

$ git clone https://github.com/mkatouda/pystrainfilter.git
$ cd pystrainfilter
$ pip install .
```

## Command line usage help
```
$ pystrainfilter -h
sage: pystrainfilter [-h] [-i INPUT] [-s SDF_PATH] [--emax-total-strain EMAX_TOTAL_STRAIN] [--emax-torsion EMAX_TORSION] [--scriptpath SCRIPT_PATH] [-v]

StrainFilter system call interface

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        yaml style input file
  -s SDF_PATH, --sdf SDF_PATH
                        input sdf file
  --emax-total-strain EMAX_TOTAL_STRAIN
                        maximum total strain energy
  --emax-torsion EMAX_TORSION
                        maximum dihedral torsion strain energy
  --scriptpath SCRIPT_PATH
                        install path of StrainFilter
  -v, --verbose         Verbose option

```

## Run sample job
### Edit bash script(sf.sh): replace scriptpath value to pato to your StrainFilter installed directory
```
$ cd /path/to/pystrainfilter
$ cp –rp sample test && cd test
$ ls
lig0_vinaout.sdf  sf.sh
$vi sf.sh
#!/bin/bash

. ${HOME}/miniconda3/etc/profile.d/conda.sh
conda activate py310-strainfilter

sdf='lig0_vinaout.sdf'
emaxtotalstrain=6.5
emaxtorsion=3.0
# Edit scriptpath to your StrainFIlter installed dirctory
scriptpath=${HOME}/data/scripts/ChemInfTools/apps/strainfilter

pystrainfilter -s ${sdf} --emax-total-strain ${emaxtotalstrain} \
    --emax-torsion ${emaxtorsion} --scriptpath ${scriptpath}
```

### run test job
```
$ bash ./sf.sh

```