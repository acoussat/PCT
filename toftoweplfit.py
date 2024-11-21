import json
import sys

import numpy as np

import itk
from itk import RTK as rtk

def poly2(a, b, c, d, e, f):
    return lambda x, y: a * x**2 + b * x * y + c * y**2 + d * x + e * y + f

def make_hull(quadric, angle):
    hull = rtk.QuadricShape.New()
    hull.SetA(quadric[0])
    hull.SetB(quadric[1])
    hull.SetC(quadric[2])
    hull.SetD(quadric[3])
    hull.SetE(quadric[4])
    hull.SetF(quadric[5])
    hull.SetG(quadric[6])
    hull.SetH(quadric[7])
    hull.SetI(quadric[8])
    hull.SetJ(quadric[9])

    theta = np.deg2rad(angle)
    hull.Rotate(itk.matrix_from_array([
        [np.cos(theta), 0., np.sin(theta)],
        [0., 1., 0.],
        [-np.sin(theta), 0., np.cos(theta)]
    ]))

    return hull

def tof_to_wepl(fit, quadric, angle):

    hull = make_hull(quadric, angle)

    try:
        with open(fit, 'r', encoding='utf-8') as f:
            coeffs_2d = json.load(f)
            poly2d = poly2(*coeffs_2d)

            def _wepl_to_tof(p_u, d_u, p_d, d_d, tof):

                # Air distance between upstream detector and object: d_uo, TOF: tof_uo
                # In the phantom: d_o, tof_o
                # Between the object and downstream detector: d_od, tof_od
                # Measured TOF = tof_uo+tof_o+tof_od
                # tof_uo is a function of energy (known) and distance (known), so it is known
                # To compute tof_uo, we first compute d_uo, then put it in the 2D fit with WEPL=0
                # For a given d_uo and a WEPL=0, we get a 1D polynomial for which the roots correspond to TOF between upstream detector and object
                # Therefore we can infer tof_o+tof_od as TOF-tof_uo
                # From there, we can use the fit with d_od (known) and tof_o+tof_od (known)
                # This gives us the WEPL

                d_d_r = [-d for d in d_d]

                intersect_u, d_uo, _ = hull.IsIntersectedByRay(p_u, d_u)  # assuming d_u is unitary
                intersect_d, d_od, _ = hull.IsIntersectedByRay(p_d, d_d_r)
                if not (intersect_u and intersect_d):
                    return 0.

                coeffs_1d = [coeffs_2d[2], coeffs_2d[1] * d_uo + coeffs_2d[4], coeffs_2d[0] * d_uo**2 + coeffs_2d[3] * d_uo + coeffs_2d[5]]  # y=[c,bx+e,ax2+dx+f]
                roots = np.roots(coeffs_1d)
                tof_uo = roots[1]  # check if always the correct one?

                tof_o_plus_tof_od = tof - tof_uo
                return poly2d(d_od, tof_o_plus_tof_od)

            return _wepl_to_tof

    except FileNotFoundError:
        sys.exit("Fit coefficient file", fit, "not found.")
