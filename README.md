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
$ conda install -c conda-forge rdkit openbabel
$ git clone https://github.com/docking-org/ChemInfTools.git
$ git clone https://github.com/mkatouda/pystrainfilter.git
$ cd pystrainfilter
$ pip install .
```

## Command line usage help
```
$ pystrainfilter -h
usage: pystrainfilter [-h] [-c COORD_PATH] [--emax-total-strain EMAX_TOTAL_STRAIN] [--emax-torsion EMAX_TORSION] [--scriptpath SCRIPT_PATH] [-v]

StrainFilter system call interface

options:
  -h, --help            show this help message and exit
  -c COORD_PATH, --coord COORD_PATH
                        input coordinate file: file format can be treated with openbabel (default: input.sdf)
  --emax-total-strain EMAX_TOTAL_STRAIN
                        maximum total strain energy (default: 6.5)
  --emax-torsion EMAX_TORSION
                        maximum dihedral torsion strain energy (default: 1.8)
  --scriptpath SCRIPT_PATH
                        install path of StrainFilter (default: ./)
  -v, --verbose         Verbose option (default: False)
```

## Run sample job
### Edit bash script(sf.sh): replace scriptpath value to pato to your StrainFilter installed directory
```
$ cd /path/to/pystrainfilter
$ cp -rp sample test && cd test
$ ls
lig0_vinaout.sdf  sf.sh
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
    --emax-torsion ${emaxtorsion} --scriptpath ${scriptpath}
```

### run test job
```
$ bash ./sf.sh
```

## Run job with multiple input files
```
#!/bin/bash                                                                                                                                                                     

. ${HOME}/miniconda3/etc/profile.d/conda.sh
conda activate py310-strainfilter

inext='.pdbqt'
emaxtotalstrain=6.5
emaxtorsion=3.0
# Edit scriptpath to your StrainFIlter installed dirctory                                                                                                                       
scriptpath=${HOME}/data/scripts/ChemInfTools/apps/strainfilter

for coord in *${inext}; do
    echo ${coord}
    pystrainfilter -c ${coord} --emax-total-strain ${emaxtotalstrain} \
        --emax-torsion ${emaxtorsion} --scriptpath ${scriptpath}
done
```
