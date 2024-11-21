#!/usr/bin/env python
import argparse
import os
import json
from multiprocessing import Pool

import matplotlib.pyplot as plt
import opengate as gate
import uproot
import numpy as np
from scipy.optimize import curve_fit

DEFAULT_OUTPUT = 'pcttoffit'
DEFAULT_NUMBER_OF_PARTICLES = int(1e4)
DEFAULT_DETECTOR_SPACING = 220.  # mm
DEFAULT_AIR_SAMPLES = int(DEFAULT_DETECTOR_SPACING // 10)
DEFAULT_WATER_SAMPLES = DEFAULT_AIR_SAMPLES


def poly2(a, b, c, d, e, f):
    return lambda x, y: a * x**2 + b * x * y + c * y**2 + d * x + e * y + f


def poly2_fit(xy, a, b, c, d, e, f):
    x = xy[:, 0]
    y = xy[:, 1]
    return poly2(a, b, c, d, e, f)(x, y)


def get_air_lengths(detector_spacing, air_samples):
    return np.linspace(0., detector_spacing, air_samples)


def get_water_lengths(air_length, detector_spacing, water_samples):
    max_water_length = detector_spacing - air_length
    return np.linspace(0., max_water_length, water_samples)


def output_format(air_length, water_length, output=DEFAULT_OUTPUT):
    return os.path.join(output, f'air_{int(air_length)}', f'water_{int(water_length)}')


def tof_fit2_mc(
    water_length_mm,
    air_length_mm,
    output=DEFAULT_OUTPUT,
    phantom_width_cm=40,
    number_of_particles=1e4,
    visu=False,
    verbose=False
):

    if verbose:
        print(f"Running MC for {air_length_mm} mm of air and {water_length_mm} mm of water…")
        print(str(locals()))

    source_energy_mev = 200

    output_base = output_format(air_length_mm, water_length_mm, output)

    out_position_mm = 0.  # arbitrary
    in_position_mm = out_position_mm - air_length_mm - water_length_mm

    # Units
    nm = gate.g4_units.nm
    mm = gate.g4_units.mm
    cm = gate.g4_units.cm
    m = gate.g4_units.m
    sec = gate.g4_units.second
    MeV = gate.g4_units.MeV

    # Simulation
    sim = gate.Simulation()

    sim.random_engine = 'MersenneTwister'
    sim.random_seed = 'auto'
    sim.run_timing_intervals = [[0 * sec, 1 * sec]]
    sim.check_volumes_overlap = False
    sim.visu = visu
    sim.visu_type = 'vrml'
    sim.g4_verbose = False
    sim.progress_bar = verbose
    sim.number_of_threads = 1

    # Misc
    yellow = [1, 1, 0, 1]
    blue = [0, 0, 1, 1]

    # Geometry
    sim.volume_manager.add_material_database('GateMaterials.db')
    sim.world.material = 'Air'
    sim.world.size = [1 * m, 1 * m, 1 * m]
    sim.world.color = [0, 0, 0, 0]

    # Phantom
    if water_length_mm > 0.:
        phantom = sim.add_volume('Box', name='Phantom')
        phantom.size = [phantom_width_cm * cm, phantom_width_cm * cm, water_length_mm * mm]
        phantom.translation = [0. * cm, 0. * cm, (-air_length_mm - water_length_mm / 2) * mm]
        phantom.material = 'Water'
        phantom.color = blue

    # Beam
    source = sim.add_source('GenericSource', 'mybeam')
    source.particle = 'proton'
    source.energy.mono = source_energy_mev * MeV
    source.energy.type = 'mono'
    source.position.type = 'box'
    source.position.size = [1 * nm, 1 * nm, 1 * nm]
    source.position.translation = [0 * mm, 0 * mm, (in_position_mm - 1.) * mm]
    source.direction.type = 'momentum'
    source.direction.momentum = [0, 0, 1]
    source.n = number_of_particles

    # Physics list
    sim.physics_manager.physics_list_name = 'G4EmStandardPhysics_option4'

    # Phase spaces

    def add_detector(name, translation, attach_to_phantom=False):
        plane = sim.add_volume('Box', 'PlanePhaseSpace' + name)
        plane.size = [phantom_width_cm * cm, phantom_width_cm * cm, 1 * nm]
        plane.translation = translation
        plane.material = 'Air'
        plane.color = yellow
        if attach_to_phantom:
            plane.mother = phantom.name

        phase_space = sim.add_actor('PhaseSpaceActor', 'PhaseSpace' + name)
        phase_space.attached_to = plane.name
        phase_space.output_filename = os.path.join(output_base, f'ps{name}.root')
        phase_space.attributes = [
            'EventID',
            'TrackID',
            'Position',
            'PreGlobalTime',
            'KineticEnergy'
        ]
        particle_filter = sim.add_filter('ParticleFilter', 'Filter' + name)
        particle_filter.particle = 'proton'

        phase_space.filters.append(particle_filter)

    add_detector('In', [0 * mm, 0 * mm, in_position_mm * mm])
    add_detector('Out', [0 * mm, 0 * mm, out_position_mm * mm])

    # Particle stats
    stat = sim.add_actor('SimulationStatisticsActor', 'stat')
    stat.output_filename = os.path.join(output_base, 'stats.txt')

    sim.run()


def tof2_fit(
    detector_spacing,
    air_samples,
    water_samples,
    output=DEFAULT_OUTPUT,
    verbose=False
):
    def print_verbose(*args, **kwargs):
        if verbose:
            print(*args, **kwargs)

    tofs = {}
    elosses = {}

    air_lengths = get_air_lengths(detector_spacing, air_samples)

    for air_length in air_lengths:

        water_lengths = get_water_lengths(air_length, detector_spacing, water_samples)

        tofs[air_length] = {}
        elosses[air_length] = {}

        for water_length in water_lengths:

            print_verbose("Processing events for air length", air_length, "mm and water length", water_length)
            output_base = output_format(air_length, water_length, output)

            tofs[air_length][water_length] = []
            elosses[air_length][water_length] = []

            try:

                data = uproot.concatenate(os.path.join(output_base, 'ps*.root'), library='np')

                if len(data['Position_Z']) == 0:
                    print_verbose("No events, skipping")
                    continue

                # Sort data if needed
                ws = data['Position_Z']
                if not np.all(ws[1:] > ws[:-1]):
                    index_sorted = np.argsort(ws)
                    for key in data.keys():
                        data[key] = data[key][index_sorted]

                for n in np.unique(data['EventID']):
                    event_mask = data['EventID'] == n

                    number_of_hits = np.sum(event_mask)
                    if number_of_hits < 2:
                        continue

                    times = data['PreGlobalTime'][event_mask]
                    tof = times[1] - times[0]

                    eloss = data['KineticEnergy'][event_mask][0] - data['KineticEnergy'][event_mask][1]

                    tofs[air_length][water_length].append(tof)
                    elosses[air_length][water_length].append(eloss)

            except FileNotFoundError:
                continue

    def fit(ys, ylabel, filename):

        def points_gen():
            for air_length in air_lengths:
                for water_length in get_water_lengths(air_length, detector_spacing, water_samples):
                    median = np.median(ys[air_length][water_length])
                    if np.isnan(median):
                        continue
                    yield [air_length, median, water_length]

        points = np.array(list(points_gen()))
        x_air = points[:, 0]
        y = points[:, 1]
        z_wepl = points[:, 2]

        xy = list(zip(x_air, y))
        p, *_ = curve_fit(poly2_fit, xy, z_wepl)
        print_verbose("Fitted parameters:", p)

        with open(os.path.join(output, filename + '.dat'), 'w', encoding='utf-8') as f:
            json.dump(list(p), f)
        with open(os.path.join(output, filename + '_points.dat'), 'w', encoding='utf-8') as f:
            for xx, yy, zz in zip(x_air, y, z_wepl):
                f.write(f'{xx} {yy} {zz}\n')

        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

        ax.scatter(x_air, y, z_wepl, marker='+')
        ax.set_xlabel("Air length [mm]")
        ax.set_ylabel(ylabel)
        ax.set_zlabel("WEPL [mm]")

        X = np.linspace(x_air.min(), x_air.max(), 100)
        Y = np.linspace(y.min(), y.max(), 100)
        X, Y = np.meshgrid(X, Y)
        Z = poly2(*p)(X, Y)
        ax.plot_surface(X, Y, Z, color='gray', alpha=.5)

        plt.savefig(os.path.join(output, filename + '.pdf'))
        plt.show()

    fit(elosses, "Energy loss [MeV]", 'eloss_fit')
    fit(tofs, "TOF [ns]", 'tof_fit')


def pcttoffit(
    output,
    number_of_particles,
    detector_spacing,
    air_samples,
    water_samples,
    phantom_width,
    visu,
    verbose
):

    with Pool(maxtasksperchild=1) as pool:

        results = []

        air_lengths = get_air_lengths(detector_spacing, air_samples)
        for air_length in air_lengths:

            if air_length < detector_spacing:
                water_lengths = get_water_lengths(air_length, detector_spacing, water_samples)
            else:
                water_lengths = [0.]

            for water_length in water_lengths:

                if air_length == 0. and water_length == 0.:
                    continue

                result = pool.apply_async(tof_fit2_mc, kwds={
                    'water_length_mm': water_length,
                    'air_length_mm': air_length,
                    'output': output,
                    'phantom_width_cm': phantom_width,
                    'number_of_particles': number_of_particles,
                    'visu': visu,
                    'verbose': verbose
                })
                results.append(result)

        pool.close()
        pool.join()

        for result in results:
            result.wait()
            if not result.successful():
                print("Failure in MC simulation")
                exit(1)

    tof2_fit(
        detector_spacing=detector_spacing,
        air_samples=air_samples,
        water_samples=water_samples,
        output=output,
        verbose=verbose
    )


def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--output', help="Path of outputs", default=DEFAULT_OUTPUT)
    parser.add_argument('-s', '--detector-spacing', help="Detector spacing", default=DEFAULT_DETECTOR_SPACING, type=int)
    parser.add_argument('-a', '--air-samples', help="Number of air samples", default=DEFAULT_AIR_SAMPLES, type=int)
    parser.add_argument('-w', '--water-samples', help="Number of water samples", default=DEFAULT_WATER_SAMPLES, type=int)
    parser.add_argument('-n', '--number-of-particles', help="Number of generated particles", default=DEFAULT_NUMBER_OF_PARTICLES, type=int)
    parser.add_argument('--phantom-width', help="Phantom width", default=40, type=float)
    parser.add_argument('--visu', help="Visualize Monte Carlo simulation", default=False, action='store_true')
    parser.add_argument('--verbose', '-v', help="Verbose execution", default=False, action='store_true')
    args_info = parser.parse_args()

    pcttoffit(**vars(args_info))

if __name__ == '__main__':
    main()
