#!/usr/bin/env python
import argparse
import uproot
import itk

from collections import OrderedDict

import numpy as np
import numpy.lib.recfunctions as rfn

def join_phase_spaces(ps_in, ps_out, no_nuclear):

    merge_columns = ['RunID', 'EventID']
    if no_nuclear:
        merge_columns.append('TrackID')

    ps_in_unique_columns = np.delete(ps_in.dtype.names, np.argwhere(np.isin(ps_in.dtype.names, merge_columns)))
    ps_out_unique_columns = np.delete(ps_out.dtype.names, np.argwhere(np.isin(ps_out.dtype.names, merge_columns)))

    pairs_names = [
        *merge_columns,
        *[t + '_in' for t in ps_in_unique_columns],
        *[t + '_out' for t in ps_out_unique_columns]
    ]

    pairs_types = [
        *[ps_in.dtype[t] for t in merge_columns],
        *[ps_in.dtype[t] for t in ps_in_unique_columns],
        *[ps_out.dtype[t] for t in ps_out_unique_columns]
    ]

    pairs = []
    ps_in_it = iter(dict(zip(ps_in.dtype.names, x)) for x in ps_in)
    ps_out_it = iter(dict(zip(ps_out.dtype.names, x)) for x in ps_out)

    ps_in_elem = next(ps_in_it)
    ps_out_elem = next(ps_out_it)

    while True:
        try:

            while ps_in_elem['RunID'] < ps_out_elem['RunID']:
                ps_in_elem = next(ps_in_it)

            while ps_in_elem['EventID'] < ps_out_elem['EventID'] and ps_in_elem['RunID'] == ps_out_elem['RunID']:
                ps_in_elem = next(ps_in_it)
            
            if no_nuclear:
                while ps_in_elem['TrackID'] < ps_out_elem['TrackID'] and ps_in_elem['EventID'] == ps_out_elem['EventID'] and ps_in_elem['RunID'] == ps_out_elem['RunID']:
                    ps_in_elem = next(ps_in_it)

            if (ps_in_elem['RunID'] == ps_out_elem['RunID']) and (ps_in_elem['EventID'] == ps_out_elem['EventID']) and ((not no_nuclear) or ps_in_elem['TrackID'] == ps_out_elem['TrackID']):
                pair = (
                    *[ps_in_elem[k] for k in merge_columns],
                    *[ps_in_elem[k] for k in ps_in_unique_columns],
                    *[ps_out_elem[k] for k in ps_out_unique_columns]
                )
                pairs.append(pair)

            ps_out_elem = next(ps_out_it)

        except StopIteration:
            break

    pairs_dtypes = list(zip(pairs_names, pairs_types))
    return np.array(pairs, dtype=pairs_dtypes)


def main():

    parser = argparse.ArgumentParser(description="Pair corresponding protons from GATE ROOT files")
    parser.add_argument('-i', '--input-in', help="Root phase space file of particles before object", required=True)
    parser.add_argument('-j', '--input-out', help="Root phase space file of particles after object", required=True)
    parser.add_argument('-o', '--output', help="Output file name", required=True)
    parser.add_argument('--plane-in', help="Plane position of incoming protons", required=True, type=float)
    parser.add_argument('--plane-out', help="Plane position of outgoing protons", required=True, type=float)
    parser.add_argument('--min-run', help="Minimum run (inclusive)", default=0, type=int)
    parser.add_argument('--max-run', help="Maximum run (exclusive)", default=1e6, type=int)
    parser.add_argument('--no-nuclear', help="Remove inelastic nuclear collisions", default=False, action='store_true')
    parser.add_argument('--verbose', '-v', help="Verbose execution", default=False, action='store_true')
    parser.add_argument('--proju', help="Provide the name of the first axis in the root file", default='Y')
    parser.add_argument('--projv', help="Provide the name of the second axis in the root file", default='Z')
    parser.add_argument('--projw', help="Provide the name of the third axis in the root file", default='X')
    parser.add_argument('--wweight', help="Weight of the third axis", default=-1.)
    parser.add_argument('--psin', help="Name of tree in input phase space", default='PhaseSpace')
    parser.add_argument('--psout', help="Name of tree in output phase space", default='PhaseSpace')
    args_info = parser.parse_args()

    if args_info.verbose:
        def verbose(message):
            print(message)
    else:
        def verbose(message):
            pass

    def load_tree_as_df(root_file, tree_name):

        branch_names = [
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
        ]

        tree = uproot.open(root_file)[tree_name]
        branches = tree.arrays(branch_names, library='np')

        dtype = [(name, branch.dtype) for name, branch in branches.items()]
        ps = np.rec.recarray((len(branches['RunID']), ), dtype=dtype)
        for branch_name in branch_names:
            ps[branch_name] = branches[branch_name]

        ps = rfn.rename_fields(ps, {
            'Position_' + str(args_info.proju): 'u',
            'Position_' + str(args_info.projv): 'v',
            'Position_' + str(args_info.projw): 'w',
        })
        ps = rfn.rename_fields(ps, {
            'Direction_' + str(args_info.proju): 'du',
            'Direction_' + str(args_info.projv): 'dv',
            'Direction_' + str(args_info.projw): 'dw',
        })

        ps = ps[(ps['RunID'] >= args_info.min_run) & (ps['RunID'] < args_info.max_run)]

        return ps

    ps_in = load_tree_as_df(args_info.input_in, 'PhaseSpaceIn')
    ps_in['w'] = args_info.plane_in
    verbose("Read input phase space:\n" + str(ps_in))

    ps_out = load_tree_as_df(args_info.input_out, 'PhaseSpaceOut')
    ps_out['w'] = args_info.plane_out
    verbose("Read output phase space:\n" + str(ps_out))

    merge_columns = ['RunID', 'EventID']
    if args_info.no_nuclear:
        merge_columns.append('TrackID')
    
    # Remove duplicates
    _, unique_index = np.unique(ps_in[merge_columns], return_index=True)
    ps_in = ps_in[unique_index]

    # The code below does not work because join_by does not handle duplicates, but we want to be able to match several output protons with one single input proton
    # pairs = rfn.join_by(
    #     merge_columns,
    #     ps_in,
    #     ps_out,
    #     jointype='inner',
    #     r1postfix='_in',
    #     r2postfix='_out',
    #     usemask=False,
    #     asrecarray=True
    # )
    # Instead, we must rely on some "hand-crafted" join_by function
    pairs = join_phase_spaces(ps_in, ps_out, args_info.no_nuclear)

    verbose("Merged input and output phase spaces:\n" + str(pairs))

    number_of_runs = pairs['RunID'].max() + 1
    verbose("Identified number of runs: " + str(number_of_runs))

    ComponentType = itk.ctype('float')
    PixelType = itk.Vector[ComponentType, 3]
    ImageType = itk.Image[PixelType, 2]

    run_range = range(args_info.min_run, min(number_of_runs, args_info.max_run))
    for r in run_range:
        ps_run = pairs[pairs['RunID'] == r]
        if len(df_run) == 0:
            continue

        ps_np = np.empty(shape=(len(ps_run), 5, 3), dtype=np.float32)
        ps_np[:,0,0] = ps_run['u_in']
        ps_np[:,0,1] = ps_run['v_in']
        ps_np[:,0,2] = ps_run['w_in']
        ps_np[:,1,0] = ps_run['u_out']
        ps_np[:,1,1] = ps_run['v_out']
        ps_np[:,1,2] = ps_run['w_out']
        ps_np[:,2,0] = ps_run['du_in']
        ps_np[:,2,1] = ps_run['dv_in']
        ps_np[:,2,2] = ps_run['dw_in'] * args_info.wweight
        ps_np[:,3,0] = ps_run['du_out']
        ps_np[:,3,1] = ps_run['dv_out']
        ps_np[:,3,2] = ps_run['dw_out'] * args_info.wweight
        ps_np[:,4,0] = ps_run['KineticEnergy_in']
        ps_np[:,4,1] = ps_run['KineticEnergy_out']
        ps_np[:,4,2] = ps_run['TrackID'] if args_info.no_nuclear else ps_run['TrackID_out']

        df_itk = itk.GetImageFromArray(ps_np, ttype=ImageType)

        output_file = args_info.output.replace('.', f'{r:04d}.')
        itk.imwrite(df_itk, output_file)
        verbose(f"Wrote file {output_file}.")

if __name__ == '__main__':
  main()