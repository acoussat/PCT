#!/usr/bin/env python
import argparse
import json
import sys
import numpy as np
import itk
from itk import PCT as pct
from itk import RTK as rtk
from toftoweplfit import tof_to_wepl, make_hull

m0 = 938.27208943  # MeV/c2
c = 299.792458  # mm/ns

LMComponentType = itk.ctype('float')
LMPixelType = itk.Vector[LMComponentType, 3]
LMImageType = itk.Image[LMPixelType, 2]

RSPPixelType = itk.ctype('double')
RSPImageType = itk.Image[RSPPixelType, 3]

def rsp_interp(rsp_itk):

    origin = rsp_itk.GetOrigin()
    index = rsp_itk.GetLargestPossibleRegion().GetIndex()
    size = rsp_itk.GetLargestPossibleRegion().GetSize()
    spacing = rsp_itk.GetSpacing()
    liif = itk.LinearInterpolateImageFunction.New(rsp_itk)

    lower_bounds = [index[i] + origin[i] - spacing[i] / 2 for i in range(3)]
    upper_bounds = [index[i] + origin[i] + size[i] + spacing[i] / 2 for i in range(3)]

    def _rsp_interp(p):
        # We check manually if cindex is within the image bounds, because wrapping seems broken
        if all(lower_bounds[i] <= p[i] <= upper_bounds[i] for i in range(3)):
            rsp = liif.Evaluate(p)
            if rsp < 0.:
                # Can happen when processing a reconstructed RSP
                return 0.
            return rsp
        return 0.

    return _rsp_interp

def velocity(e):
    # Same equation as Ulrich-Pur except that c is not in natural units
    # The terms c² cancel out
    return c * (e / (e + m0)) * np.sqrt(1 + 2 * (m0 / e))

def tof_in_vacuum(p1, p2, e):
    u1, v1, w1 = p1
    u2, v2, w2 = p2

    d = np.sqrt((u1 - u2)**2 + (v1 - v2)**2 + (w1 - w2)**2)
    v = velocity(e)  # e is constant in vacuum

    return d / v

def compute_tof(pos_in, pos_out, dir_in, dir_out, sp_fit, rsp, polydeg, e0, dw, quadric=None, angle=None, skip_vacuum=False):

    tof_tot = 0.

    pos_in_q = pos_in
    pos_out_q = pos_out
    w_out = pos_out[2]
    if skip_vacuum:
        hull = make_hull(quadric, angle)
        dir_out_r = [-d for d in dir_out]
        intersect_in, d_in, _ = hull.IsIntersectedByRay(pos_in, dir_in)
        intersect_out, d_out, _ = hull.IsIntersectedByRay(pos_out, dir_out_r)
        if intersect_in != intersect_out:
            # The path intersects the object on only one end?
            # This would not make sense
            sys.exit("The quadric seems malformed, exiting")
        if intersect_in:
            pos_in_q = [pos_in[i] + d_in * dir_in[i] for i in range(3)]
            tof_tot += tof_in_vacuum(pos_in, pos_in_q, e0)
        if intersect_out:
            pos_out_q = [pos_out[i] + d_out * dir_out_r[i] for i in range(3)]
            w_out = pos_out_q[2]

        if not intersect_in and not intersect_out:
            # If there is not intersection, the proton will never intersect with the object
            # But we still have to compute the TOF in air
            return tof_in_vacuum(pos_in, pos_out, e0)

    # Loop variables
    e = e0
    w = pos_in_q[2]

    mlp = pct.PolynomialMLPFunction.New()
    mlp.SetPolynomialDegree(polydeg)
    mlp.Init(pos_in_q, pos_out_q, dir_in, dir_out)

    while w < w_out:

        u, v = mlp.Evaluate([w])
        u, v = u[0], v[0]

        w2 = w + dw
        u2, v2 = mlp.Evaluate([w2])
        u2, v2 = u2[0], v2[0]

        d = np.sqrt((u - u2)**2 + (v - v2)**2 + (w - w2)**2)

        rsp_p = rsp([u, v, w])
        rsp_p2 = rsp([u2, v2, w2])
        rsp_avg = (rsp_p + rsp_p2) / 2
        sp_h2o = np.polyval(sp_fit, e)  # From fit
        sp = rsp_avg * sp_h2o
        eloss = sp * d
        e2 = e - eloss

        if e2 < 0.:
            # The proton has lost all its energy and will never reach the downstream detector
            return np.inf

        e_avg = (e + e2) / 2
        v = velocity(e_avg)

        tof = d / v
        tof_tot += tof

        e = e2
        w = w2

    if skip_vacuum and w < pos_out[2]:
        # if we stopped before the exit detector, then we stopped at the object ending
        # in that case, we compute the TOF in the vacuum after the object
        u, v = mlp.Evaluate([w])
        u, v = u[0], v[0]
        tof_tot += tof_in_vacuum(pos_out, [u, v, w], e)

    return tof_tot

def pcttofforward(
    input_listmode,
    output_listmode,
    sp_fit,
    rsp_map,
    initial_energy,
    step,
    tof_fit,
    quadric,
    angle,
    polydeg,
    skip_vacuum=False,
    verbose=False
):
    if verbose:
        print(f"Reading input file {input_listmode}…")
    lm_itk = itk.imread(input_listmode)
    lm_np = itk.GetArrayFromImage(lm_itk)

    ttw = tof_to_wepl(tof_fit, quadric, angle)

    rsp_itk = itk.imread(rsp_map)
    rsp = rsp_interp(rsp_itk)

    with open(sp_fit, 'r', encoding='utf-8') as file:
        sp_fit = json.load(file)

    for i, lm_row in enumerate(lm_np):
        pos_in = [float(v) for v in lm_row[0, :]]
        pos_out = [float(v) for v in lm_row[1, :]]
        dir_in = [float(v) for v in lm_row[2, :]]
        dir_out = [float(v) for v in lm_row[3, :]]

        tof = compute_tof(pos_in, pos_out, dir_in, dir_out, sp_fit, rsp, polydeg, initial_energy, step, quadric, angle, skip_vacuum)
        if tof == np.inf:
            # The proton has lost all its energy, so we ignore it
            wepl = 0.
        else:
            # Convert TOF to wepl
            wepl = ttw(pos_in, dir_in, pos_out, dir_out, tof)

        lm_row[4, 0] = 0.
        lm_row[4, 1] = wepl

        if verbose and i % 1000 == 0:
            print(f"{i} protons processed")

    if verbose:
        print(f"Writing output file {output_listmode}…")
    lm_out_itk = itk.GetImageFromArray(lm_np, ttype=LMImageType)
    itk.imwrite(lm_out_itk, output_listmode)

def main():

    parser = argparse.ArgumentParser(description="Apply TOF forward model to PCT list-mode data")
    parser.add_argument('-i', '--input-listmode', help="Input list-mode data", required=True)
    parser.add_argument('-o', '--output-listmode', help="Output list-mode data", required=True)
    parser.add_argument('--sp-fit', help="Fit of stopping power of water", required=True)
    parser.add_argument('-r', '--rsp-map', help="Relative stopping power map", required=True)
    parser.add_argument('-e', '--initial-energy', help="Initial energy of the protons in MeV", default=200.)
    parser.add_argument('--step', help="Step size (in mm along z)", type=float, default=1.)
    parser.add_argument('--tof-fit', help="Fit for TOF calibration")
    parser.add_argument('--quadric', help="Quadric representing the hull of the object", type=float, nargs=10, required=False)
    parser.add_argument('--angle', help="Angle of the projection", type=float, required=True)
    parser.add_argument('--polydeg', help="Polynomial degree for MLP estimation", type=int, default=5)
    parser.add_argument('--skip-vacuum', help="Skip vacuum when computing TOF (should be equivalent, about 30%% faster)", action='store_true')
    parser.add_argument('--verbose', '-v', help="Verbose execution", default=False, action='store_true')
    args_info = parser.parse_args()

    pcttofforward(**vars(args_info))

if __name__ == '__main__':
    main()