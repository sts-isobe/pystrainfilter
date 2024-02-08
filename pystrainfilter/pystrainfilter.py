#!/usr/bin/env python

import os
import glob
import argparse
import subprocess

from openbabel import pybel
from rdkit import Chem
import numpy as np
import pandas as pd
try:
    import yaml
except:
    pass


scname = ['total_strain', 'dihedral_torsion_strain']

def get_parser():

    class customHelpFormatter(argparse.ArgumentDefaultsHelpFormatter,
                              argparse.RawTextHelpFormatter):
        pass

    parser = argparse.ArgumentParser(
        description="StrainFilter system call python interface",
        formatter_class=customHelpFormatter
    )
    parser.add_argument(
        "-i","--inp",
        help="yaml style input file",
        type=str,
        default=None
    )
    parser.add_argument(
        "-c","--coord",
        help="input coordinate file or directory: file format can be treated with openbabel",
        type=str,
        dest='coord',
        default='input.sdf'
    )
    parser.add_argument(
        "--emax-total-strain",
        help="maximum total strain energy",
        type=float,
        dest='emax_total_strain',
        default=6.0
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
        default='./ChemInfTools/apps/strainfilter'
    )
    parser.add_argument(
        "-o","--out",
        help="basename of output csv file",
        type=str,
        dest='out',
        default='strain'
    )
    parser.add_argument(
        "-s","--savesdf",
        help="save sdf file of filtered structures",
        action='store_true',
        dest='savesdf'
    )
    parser.add_argument(
        "--addH",
        help="add hydrogen atoms to filtered structures in sdf file",
        action='store_true',
        dest='addH'
    )
    parser.add_argument(
        "--savescr",
        help="save scracth files option",
        action='store_true',
        dest='savescr'
    )
    parser.add_argument(
        "-v","--verbose",
        help="Verbose option",
        action='store_true',
        dest='verbose'
    )
    return parser.parse_args()


def set_config(args):

    # Read config yaml file
    if args.inp is not None and os.path.isfile(args.inp):
        with open(args.inp, 'r') as f:
            conf = yaml.safe_load(f)
    else:
        conf = {}

    # Set up default config values from program arguments
    conf_def = vars(args).copy()
    [conf.setdefault(k, v) for k, v in conf_def.items()]

    return conf


def run_strainfilter(coord_path, script_path, emax_total_strain=6.5,
                     emax_torsion=1.8, savescr=False, verbose=False, **kwargs):

    basename, ext = os.path.splitext(os.path.basename(coord_path))
    print('input coordinate file:', coord_path, flush=True)
    ext = ext.lower()
    mol2_path = basename + '.mol2'
    #ob_mols = pybel.readfile(ext[1:], coord_path)
    #mol2file = pybel.Outputfile('mol2', mol2_path, overwrite=True)
    #for ob_mol in ob_mols:
    #    mol2file.write(ob_mol)
    #mol2file.close()
    cmd = [
        'obabel', '-i'+ext[1:], coord_path, '-omol2', '-O', mol2_path, '-xu'
    ]
    print(
        'Convert coordinate file using openbabel:\n' + ' '.join(cmd),
        flush=True
    )
    results = subprocess.run(
        cmd, capture_output=True, check=True, text=True
    )
    if verbose: print(results.stdout, results.stderr, flush=True)

    cw_dir = os.getcwd()
    os.chdir(script_path)

    cmd = [
        'python', 'Torsion_Strain.py', cw_dir + '/' + mol2_path
    ]
    print(
        'Run StrainFilter:\n' + ' '.join(cmd), flush=True
    )

    results = subprocess.run(
        cmd, capture_output=True, check=True, text=True
    )
    if verbose: print(results.stdout, results.stderr, flush=True)

    os.chdir(cw_dir)

    csv_path = basename + '_Torsion_Strain.csv'
    try:
        df = pd.read_csv(csv_path, header=None)[[1, 5]]
        df.index.name = 'SID'
        df.columns = scname
        print('StrainFilter output csv_path:', csv_path, '\n', df, '\n', flush=True)
        name = [os.path.basename(coord_path) for i in range(len(df))]
        score = df.to_numpy()
    except:
        name = [os.path.basename(coord_path)]
        score = np.zeros([1, 2])
        score[:, :] = np.nan
        pass

    if not savescr:
        for file_path in [mol2_path, csv_path]:
            if os.path.exists(file_path): os.remove(file_path)

    if verbose: print('name:', name, 'score:', score, flush=True)

    return name, score


def sdf_write(df, coord_dir_path='.', inext='.sdf', out='strain', addH=False,
              verbose=False, **kwargs):

    if verbose: print('Write structures to a sdf file.', flush=True)

    sdfout_path = out + '_geom_filtered.sdf'
    with pybel.Outputfile('sdf', sdfout_path, overwrite=True) as out1:

        filenames = df['filename'].unique().tolist()
        #for idx, item in df.iterrows():
        for filename in filenames:
            infilename = '{}/{}'.format(coord_dir_path, filename)
            if verbose: print('infilename:', infilename, flush=True)
            mols = list(pybel.readfile(inext[1:], infilename))
            #m = next(pybel.readfile(inext[1:], infilename))
            df_s = df.query('filename == "{}"'.format(filename))
            for idx, item in df_s.iterrows():
                if verbose: print(idx, '\n', item, flush=True)
                m = mols[item['SID']]
                if addH:
                    with pybel.Outputfile('pdb', 'tmp.pdb', overwrite=True) as out2:
                        out2.write(m)
                        m = next(pybel.readfile('pdb', 'tmp.pdb'))
                        m.OBMol.AddHydrogens()
                        os.remove('tmp.pdb')
                m.data['input_file_name'] = item['filename']
                m.data['structure_ID'] = item['SID']
                for scn in scname:
                    m.data[scn] = item[scn]
                out1.write(m)


def run_strainfilter_batch(conf):

    coord = conf['coord']
    emax_total_strain = conf['emax_total_strain']
    emax_torsion = conf['emax_torsion']
    savesdf = conf['savesdf']
    addH = conf['addH']
    verbose = conf['verbose']

    if os.path.isdir(coord):
        coord_dir_path = coord
        for ext in ['.sdf', '.pdbqt', '.mol2', '.mol']:
            coord_files ='{}/*{}'.format(coord, ext)
            coord_files = [cf for cf in sorted(glob.glob(coord_files))]
            if len(coord_files) > 0:
                inext = ext
                break
    else:
        coord_dir_path = os.path.dirname(coord)
        if len(coord_dir_path) == 0: coord_dir_path = '.'
        coord_files = [coord]
        inext = os.path.splitext(os.path.basename(coord))[1]
    ncrdfile = len(coord_files)

    if verbose:
        print('coordinate files:\n', coord_files, flush=True)
        print('coord_dir_path:', coord_dir_path, flush=True)

    name = []; sid = []; score = np.empty([0, 2])
    for idx, coord_path in enumerate(coord_files):
        print('##### File {}/{} #####'.format(idx+1, ncrdfile), flush=True)
        na, sc = run_strainfilter(coord_path, **conf)
        name += na
        score = np.concatenate([score, sc])
        sid += [i for i in range(sc.shape[0])] 

    df = pd.DataFrame(
        {
            'filename': name, 'SID': sid,
            scname[0]: score[:, 0], scname[1]: score[:, 1]
        }
    )
    csvout_path = conf['out'] + '_score_raw.csv'
    df.to_csv(csvout_path, index=False)

    df_filtered = df[
        (df[scname[0]] <= emax_total_strain) & (df[scname[1]] <= emax_torsion)
    ]
    csvout_path = conf['out'] + '_score_filtered.csv'
    df_filtered.to_csv(csvout_path, index=False)

    if savesdf:
        sdf_write(
            df=df_filtered, coord_dir_path=coord_dir_path, inext=inext, **conf
        )

    nfilt = len(df_filtered); ntot = len(df)
    ratio_filt = float(nfilt) / float(ntot)

    print(
        '\npystrainfilter finished!: #filtered/#input={}/{} ({:.3f})\n'\
        .format(nfilt, ntot, ratio_filt), flush=True
    )


def main():
    args = get_parser()
    if args.verbose: print(args, flush=True)

    conf = set_config(args)

    s = '======= Input configulations =======\n'
    for k, v in conf.items():
        s += '{}: {}\n'.format(k, v)
    s += '====================================\n'
    print(s, flush=True)

    if os.path.isdir(conf['script_path']):
        run_strainfilter_batch(conf)


if __name__ == '__main__':
    main()
