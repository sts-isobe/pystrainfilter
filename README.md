# pystrainfilter
python interface of StrainFilter

## Licence
MIT

## Install using conda yaml environment file
```
$ git clone https://github.com/mkatouda/pystrainfilter.git
$ cd pystrainfilter
$ conda create -f py310-strainfilter
$ conda activate py310-strainfilter
$ git clone https://github.com/docking-org/ChemInfTools.git
```

## Install from conda package
```
$ conda create -n py310-strainfilter -c conda-forge python=3.10
$ conda activate py310-strainfilter
$ conda install -c conda-forge rdkit openbabel pyyaml
$ git clone https://github.com/docking-org/ChemInfTools.git
$ git clone https://github.com/mkatouda/pystrainfilter.git
$ cd pystrainfilter
$ pip install .
```

## Command line usage help
```
$ pystrainfilter -h
usage: pystrainfilter [-h] [-i INP] [-c COORD] [--emax-total-strain EMAX_TOTAL_STRAIN] [--emax-torsion EMAX_TORSION] [--scriptpath SCRIPT_PATH] [-o OUT] [-s] [--addH] [--savescr] [-v]

StrainFilter system call python interface

options:
  -h, --help            show this help message and exit
  -i INP, --inp INP     yaml style input file (default: None)
  -c COORD, --coord COORD
                        input coordinate file or directory: file format can be treated with openbabel (default: input.sdf)
  --emax-total-strain EMAX_TOTAL_STRAIN
                        maximum total strain energy (default: 6.0)
  --emax-torsion EMAX_TORSION
                        maximum dihedral torsion strain energy (default: 1.8)
  --scriptpath SCRIPT_PATH
                        install path of StrainFilter (default: ./ChemInfTools/apps/strainfilter)
  -o OUT, --out OUT     basename of output csv file (default: strain)
  -s, --savesdf         save sdf file of filtered structures (default: False)
  --addH                add hydrogen atoms to filtered structures in sdf file (default: False)
  --savescr             save scracth files option (default: False)
  -v, --verbose         Verbose option (default: False)
```

## Run sample job with single sdf file
### Copy sample job directory
```
$ cd /path/to/pystrainfilter
$ cp -rp sample/sdf testsdf && cd testsdf
$ ls
lig0_vinaout.sdf  sf.sh
```

### Edit bash script(sf.sh): replace scriptpath value to path to your StrainFilter installed directory
```
$ vi sf.sh
#!/bin/bash

. ${HOME}/miniconda3/etc/profile.d/conda.sh
conda activate py310-strainfilter

coord='lig0_vinaout.sdf'
emaxtotalstrain=6.5
emaxtorsion=3.0
# Edit scriptpath to your StrainFIlter installed dirctory
scriptpath=${HOME}/data/scripts/ChemInfTools/apps/strainfilter

pystrainfilter -c ${coord} --emax-total-strain ${emaxtotalstrain} \
    --emax-torsion ${emaxtorsion} --scriptpath ${scriptpath} -s
```

### Run sample job
```
$ bash ./sf.sh
```

## Run sample job with multiple pdbqt input files
### Copy sample job directory
```
$ cd /path/to/pystrainfilter
$ cp -rp sample/pdbqt testpdbqt && cd testpdbqt
$ ls
MCULE-1005714968.pdbqt	MCULE-1037768225.pdbqt	MCULE-1062803813.pdbqt	MCULE-1089858904.pdbqt	MCULE-1089955095.pdbqt	sf.sh
```

### Edit bash script(sf.sh): replace scriptpath value to path to your StrainFilter installed directory

```
$ vi sf.sh
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
```

### Run sample job
```
$ bash ./sf.sh
```

## Using yaml input file for sample job with single sdf input file
### Copy sample job directory
```
$ cd /path/to/pystrainfilter
$ cp -rp sample/sdf testsdf && cd testsdf
$ ls
lig0_vinaout.sdf  sf.sh
```

### Prepare yaml input file
Note: yaml keyward is replaced '-' with '_' from command options.

```
$ vi input.yml
coord : 'lig0_vinaout.sdf'
emax_total_strain : 6.0
emax_torsion : 1.8
# Edit scriptpath to your StrainFIlter installed dirctory
script_path : '~/data/scripts/ChemInfTools/apps/strainfilter'
savesdf : True
```

### Make bash script(sf.sh): replace scriptpath value to path to your StrainFilter installed directory

```
$ vi sf_yml.sh
#!/bin/bash

. ${HOME}/miniconda3/etc/profile.d/conda.sh
conda activate py310-strainfilter

pystrainfilter -i input.yml
```

### Run sample job
```
$ bash ./sf_yml.sh
```
