#!/usr/bin/env python

import sys
import os
import shutil
import argparse
import subprocess

#from openbabel import pybel
from rdkit import Chem
import pandas as pd


def get_parser():
    class customHelpFormatter(argparse.ArgumentDefaultsHelpFormatter,
                              argparse.RawTextHelpFormatter):
        pass

    parser = argparse.ArgumentParser(
        description="StrainFilter system call interface"
    )
    parser.add_argument(
        "-i","--input",
        help="yaml style input file",
        type=str,
        default='input.yml'
    )
    parser.add_argument(
        "-s","--sdf",
        help="input sdf file",
        type=str,
        dest='sdf_path',
        default='input.sdf'
    )
    parser.add_argument(
        "--emax-total-strain",
        help="maximum total strain energy",
        type=float,
        dest='emax_total_strain',
        default=6.5
    )
    parser.add_argument(
        "--emax-torsion",
        help="maximum dihedral torsion strain energy",
        type=float,
        dest='emax_torsion',
        default=1.8
    )
    parser.add_argument(
        "--scriptpath",
        help="install path of StrainFilter",
        type=str,
        dest='script_path',
        default='./'
    )
    parser.add_argument(
        "-v","--verbose",
        help="Verbose option",
        action='store_true'
    )
    return parser.parse_args()


def set_config(args):
    # Read config yaml file
    #if args.inp is not None and os.path.isfile(args.inp):
    #    with open(args.inp, 'r') as f:
    #        conf = yaml.safe_load(f)
    #else:
    #    conf = {}
    conf = {}

    # Set up default config values from program arguments
    conf_def = vars(args).copy()
    [conf.setdefault(k, v) for k, v in conf_def.items()]

    return conf


def run_strainfilter(sdf_path, script_path, emax_total_strain=6.5,
                     emax_torsion=1.8, verbose=False, **kwargs):

    basename = os.path.splitext(os.path.basename(sdf_path))[0]
    mol2_path = basename + '.mol2'
    #ob_mols = pybel.readfile('sdf', sdf_path)
    #mol2file = pybel.Outputfile('mol2', mol2_path, overwrite=True)
    #for ob_mol in ob_mols:
    #    mol2file.write(ob_mol)
    #mol2file.close()
    cmd = [
        'obabel', '-isdf', sdf_path, '-omol2', '-O', mol2_path, '-xu'
    ]
    print(' '.join(cmd), flush=True)
    results = subprocess.run(
        cmd, capture_output=True, check=True, text=True
    )
    if verbose: print(results.stdout, results.stderr, flush=True)

    cw_dir = os.getcwd()
    os.chdir(script_path)

    cmd = [
        'python',
        'Torsion_Strain.py',
        cw_dir + '/' + mol2_path
    ]
    print(' '.join(cmd), flush=True)

    results = subprocess.run(
        cmd, capture_output=True, check=True, text=True
    )
    if verbose: print(results.stdout, results.stderr, flush=True)

    os.chdir(cw_dir)

    csv_path = basename + '_Torsion_Strain.csv'
    df = pd.read_csv(csv_path, header=None)[[1, 5]]
    print('csv_path:', csv_path, '\n', df, flush=True)
    df_filtered = df[(df[1] <= emax_total_strain) & (df[5] <= emax_torsion)]
    print(df_filtered, flush=True)
    idx_filtered = df_filtered.index.to_list()
    print('idx_filtered:', idx_filtered, flush=True)

    sdfout_path = basename + '_strain_filtered.sdf' 
    rd_mols = list(Chem.SDMolSupplier(sdf_path, removeHs=False))
    with Chem.SDWriter(sdfout_path) as writer:
        for i in idx_filtered:
            writer.write(rd_mols[i])

    return idx_filtered


def main():
    args = get_parser()
    if args.verbose: print(args, flush=True)

    conf = set_config(args)

    print('======= Input configulations =======', flush=True)
    for k, v in conf.items():
        print('{}: {}'.format(k, v), flush=True)
    print('====================================', flush=True)

    if os.path.isdir(conf['script_path']):
        run_strainfilter(**conf)

if __name__ == '__main__':
    main()
