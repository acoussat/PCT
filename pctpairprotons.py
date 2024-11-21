"""Merge protons from proton CT simulation.

It is assumed that ROOT files were already filtered and only contain protons.
"""

import argparse
import uproot
import itk

import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input-in', help="Root phase space file of particles before object", required=True)
parser.add_argument('-j', '--input-out', help="Root phase space file of particles after object", required=True)
parser.add_argument('-o', '--output', help="Output file name", required=True)
parser.add_argument('--plane-in', help="Plane position of incoming protons", required=True, type=float)
parser.add_argument('--plane-out', help="Plane position of outgoing protons", required=True, type=float)
parser.add_argument('--verbose', '-v', help="Verbose execution", default=False, action='store_true')
parser.add_argument('--proju', help="Provide the name of the first axis in the root file", default='Y')
parser.add_argument('--projv', help="Provide the name of the second axis in the root file", default='Z')
parser.add_argument('--projw', help="Provide the name of the third axis in the root file", default='X')
parser.add_argument('--wweight', help="Weight of the third axis", default=-1.)
parser.add_argument('--psin', help="Name of tree in input phase space", default='PhaseSpace')
parser.add_argument('--psout', help="Name of tree in output phase space", default='PhaseSpace')
args = parser.parse_args()

if args.verbose:
    def verbose(message):
        print(message)
else:
    def verbose(message):
        pass

def load_tree_as_df(root_file, tree_name):
    tree = uproot.open(root_file)[tree_name]
    df = tree.arrays([
        'RunID',
        'EventID',
        'TrackID',
        'KineticEnergy',
        'GlobalTime',
        'Position_X',
        'Position_Y',
        'Position_Z',
        'Direction_X',
        'Direction_Y',
        'Direction_Z'
    ], library='pd')
    df.rename(columns={
        'Position_' + str(args.proju): 'u',
        'Position_' + str(args.projv): 'v',
        'Position_' + str(args.projw): 'w',
    }, inplace=True)
    df.rename(columns={
        'Direction_' + str(args.proju): 'du',
        'Direction_' + str(args.projv): 'dv',
        'Direction_' + str(args.projw): 'dw',
    }, inplace=True)

    df['dw'] *= args.wweight

    return df

df_in = load_tree_as_df(args.input_in, 'PhaseSpaceIn')
df_in['w'] = args.plane_in
verbose("Read input phase space:\n" + str(df_in))

df_out = load_tree_as_df(args.input_out, 'PhaseSpaceOut')
df_out['w'] = args.plane_out
verbose("Read output phase space:\n" + str(df_out))

merge_columns = ['RunID', 'EventID']
df_pairs = pd.merge(
    df_in.drop_duplicates(merge_columns),  # drop_duplicates to mimic the behavior of C++ pctpairprotons (see comment line 272, commit 8260d3)
    df_out,
    on=merge_columns,
    how='inner',
    suffixes=('_in', '_out')
)
verbose("Merged input and output phase spaces:\n" + str(df_pairs))

max_runs = df_pairs['RunID'].max() + 1
verbose("Identified number of runs: " + str(max_runs))

ComponentType = itk.ctype('float')
PixelType = itk.Vector[ComponentType, 3]
ImageType = itk.Image[PixelType, 2]

for r in range(max_runs):
    verbose(f"Processing run {r}…")

    df_run = df_pairs[df_pairs['RunID'] == r]

    df_np = np.empty(shape=(len(df_run), 5, 3), dtype=np.float32)
    df_np[:,0,0] = df_run['u_in']
    df_np[:,0,1] = df_run['v_in']
    df_np[:,0,2] = df_run['w_in']
    df_np[:,1,0] = df_run['u_out']
    df_np[:,1,1] = df_run['v_out']
    df_np[:,1,2] = df_run['w_out']
    df_np[:,2,0] = df_run['du_in']
    df_np[:,2,1] = df_run['dv_in']
    df_np[:,2,2] = df_run['dw_in']
    df_np[:,3,0] = df_run['du_out']
    df_np[:,3,1] = df_run['dv_out']
    df_np[:,3,2] = df_run['dw_out']
    df_np[:,4,0] = df_run['KineticEnergy_in']
    df_np[:,4,1] = df_run['KineticEnergy_out']
    df_np[:,4,2] = df_run['TrackID_out']
    verbose("Converted pairs to NumPy:\n" + str(df_np))

    df_itk = itk.GetImageFromArray(df_np, ttype=ImageType)
    verbose("Converted NumPy array to ITK image.")

    output_file = args.output.replace('.', f'{r:04d}.')
    itk.imwrite(df_itk, output_file)
    verbose(f"Wrote file {output_file}.")