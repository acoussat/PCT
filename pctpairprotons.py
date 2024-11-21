#!/usr/bin/env python
import argparse

import numpy as np
import numpy.lib.recfunctions as rfn

import uproot
import itk
from itk import RTK as rtk

from toftoweplfit import tof_to_wepl

def pctpairprotons(
    input_in,
    input_out,
    output,
    plane_in,
    plane_out,
    min_run=0,
    max_run=1e6,
    no_nuclear=False,
    verbose=False,
    proju='Y',
    projv='Z',
    projw='X',
    fit=None,
    quadric=None,
    angle=None
):
    def print_verbose(*args, **kwargs):
        if verbose:
            print(*args, **kwargs)

    def load_tree_as_df(root_files):

        branch_names = [
            'RunID',
            'EventID',
            'TrackID',
            'KineticEnergy',
            'Position_X',
            'Position_Y',
            'Position_Z',
            'Direction_X',
            'Direction_Y',
            'Direction_Z',
            'PreGlobalTime'
        ]

        branches = uproot.concatenate(root_files, branch_names, library='np')

        dtype = [(name, branch.dtype) for name, branch in branches.items()]
        ps = np.rec.recarray((len(branches['RunID']), ), dtype=dtype)
        for branch_name in branch_names:
            ps[branch_name] = branches[branch_name]

        print_verbose("Sorting data…")
        np.recarray.sort(ps, order=['RunID', 'EventID', 'TrackID'])

        ps = rfn.rename_fields(ps, {
            'Position_' + str(proju): 'u',
            'Position_' + str(projv): 'v',
            'Position_' + str(projw): 'w',
        })
        ps = rfn.rename_fields(ps, {
            'Direction_' + str(proju): 'du',
            'Direction_' + str(projv): 'dv',
            'Direction_' + str(projw): 'dw',
        })

        ps = ps[(ps['RunID'] >= min_run) & (ps['RunID'] < max_run)]

        return ps

    ps_in = load_tree_as_df(input_in)
    ps_in['w'] = plane_in
    print_verbose("Read input phase space:\n" + str(ps_in))

    ps_out = load_tree_as_df(input_out)
    ps_out['w'] = plane_out
    print_verbose("Read output phase space:\n" + str(ps_out))

    merge_columns = ['RunID', 'EventID']
    if no_nuclear:
        merge_columns.append('TrackID')

    # Remove duplicates
    _, unique_index = np.unique(ps_in[merge_columns], return_index=True)
    ps_in = ps_in[unique_index]

    ps_in.dtype.names = [n if n in merge_columns else n + '_in' for n in ps_in.dtype.names]
    ps_out.dtype.names = [n if n in merge_columns else n + '_out' for n in ps_out.dtype.names]
    ps_in_uniques = [n for n in ps_in.dtype.names if n not in merge_columns]
    ps_out_uniques = [n for n in ps_out.dtype.names if n not in merge_columns]

    if no_nuclear:
        # Easy case, there should be at most one row in ps_in and ps_out for keys ['RunID', 'EventID', 'TrackID']
        intersect_u, intersect_in, intersect_out = np.intersect1d(ps_in[merge_columns], ps_out[merge_columns], return_indices=True)
        pairs = rfn.merge_arrays((intersect_u, ps_in[ps_in_uniques][intersect_in], ps_out[ps_out_uniques][intersect_out]), asrecarray=True, flatten=True)
    else:
        # More complicated, there can be more than one row in ps_out for keys ['RunID', 'EventID'] (because of 'TrackID')
        # The solution is to repeat the computation for many TrackIDs, then merge it all together
        track_max = ps_out['TrackID_out'].max()
        print_verbose("Identified maximum number of tracks: " + str(track_max))
        pairs_list = []
        for t in range(track_max + 1):
            ps_out_t = ps_out[ps_out['TrackID_out'] == t]
            intersect_u, intersect_in, intersect_out = np.intersect1d(ps_in[merge_columns], ps_out_t[merge_columns], return_indices=True)
            pairs = rfn.merge_arrays((intersect_u, ps_in[ps_in_uniques][intersect_in], ps_out_t[ps_out_uniques][intersect_out]), asrecarray=True, flatten=True)
            if len(pairs) > 0:
                pairs_list.append(pairs)
        pairs = rfn.stack_arrays(pairs_list, asrecarray=True)
        np.recarray.sort(pairs, order=['RunID', 'EventID', 'TrackID_in', 'TrackID_out'])
    print_verbose(f"Merged input and output phase spaces into {len(pairs)} pairs.")

    if fit is not None:
        print_verbose("Processing fit…")
        names = pairs.dtype.names
        for pi, pair in enumerate(pairs):
            p = dict(zip(names, pair.tolist()))
            p_u = [p['u_in'], p['v_in'], p['w_in']]
            d_u = [p['du_in'], p['dv_in'], p['dw_in']]
            p_d = [p['u_out'], p['v_out'], p['w_out']]
            d_d = [p['du_out'], p['dv_out'], p['dw_out']]
            tof = p['PreGlobalTime_out'] - p['PreGlobalTime_in']
            ttw = tof_to_wepl(fit, quadric, angle)
            wepl = ttw(p_u, d_u, p_d, d_d, tof)
            pairs[pi]['KineticEnergy_in'] = 0.
            pairs[pi]['KineticEnergy_out'] = wepl

    number_of_runs = pairs['RunID'].max() + 1
    print_verbose("Identified number of runs: " + str(number_of_runs))

    ComponentType = itk.ctype('float')
    PixelType = itk.Vector[ComponentType, 3]
    ImageType = itk.Image[PixelType, 2]

    if output is not None:
        run_range = range(min_run, min(number_of_runs, max_run))
        for r in run_range:
            ps_run = pairs[pairs['RunID'] == r]
            if len(ps_run) == 0:
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
            ps_np[:,2,2] = ps_run['dw_in']
            ps_np[:,3,0] = ps_run['du_out']
            ps_np[:,3,1] = ps_run['dv_out']
            ps_np[:,3,2] = ps_run['dw_out']
            ps_np[:,4,0] = ps_run['KineticEnergy_in']
            ps_np[:,4,1] = ps_run['KineticEnergy_out']
            ps_np[:,4,2] = ps_run['TrackID'] if no_nuclear else ps_run['TrackID_out']

            df_itk = itk.GetImageFromArray(ps_np, ttype=ImageType)

            output_file = output.replace('.', f'{r:04d}.')
            itk.imwrite(df_itk, output_file)
            print_verbose(f"Wrote file {output_file}.")

    return pairs

def main():

    parser = argparse.ArgumentParser(description="Pair corresponding protons from GATE ROOT files")
    parser.add_argument('-i', '--input-in', help="Root phase space files of particles before object", required=True, nargs='+')
    parser.add_argument('-j', '--input-out', help="Root phase space files of particles after object", required=True, nargs='+')
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
    parser.add_argument('--fit', help="JSON file that contains data fit")
    parser.add_argument('--quadric', help="Quadric representing the hull of the object", type=float, nargs=10, required=False)
    parser.add_argument('--angle', help="Angle of the current projection", type=float, required=False)
    args_info = parser.parse_args()

    pctpairprotons(**vars(args_info))

if __name__ == '__main__':
    main()
