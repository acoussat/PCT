#!/usr/bin/env python
import argparse
import os
import json
from multiprocessing import Pool

import matplotlib.pyplot as plt
import opengate as gate
import uproot
import numpy as np

DEFAULT_OUTPUT = 'pctspwatercalibration'
DEFAULT_ENERGY_MAX = 200.  # MeV
DEFAULT_ENERGY_MIN = 10.  # MeV
DEFAULT_ENERGY_SAMPLES = 10
DEFAULT_NUMBER_OF_PARTICLES = 50000


def output_format(output, energy):
    return os.path.join(output, f'{int(energy)}MeV')


def pctwatercalibration_mc(
    energy,
    number_of_particles=DEFAULT_NUMBER_OF_PARTICLES,
    output=DEFAULT_OUTPUT,
    visu=False,
    verbose=False
):

    if verbose:
        print("Running MC for energy of", energy, "MeV")
        print(str(locals()))

    output_base = output_format(output, energy)

    # Units
    nm = gate.g4_units.nm
    mm = gate.g4_units.mm
    cm = gate.g4_units.cm
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
    blue = [0, 0, 1, .25]

    # Geometry
    sim.volume_manager.add_material_database('GateMaterials.db')
    sim.world.material = 'Air'
    sim.world.size = [1. * cm, 1. * cm, 1. * cm]
    sim.world.color = [0, 0, 0, 0]

    # Phantom
    phantom = sim.add_volume('Box', name='Phantom')
    phantom.size = [1. * mm, 1. * mm, 1. * mm]
    phantom.translation = [0. * cm, 0. * cm, 0. * cm]
    phantom.material = 'Water'
    phantom.color = blue

    # Beam
    source = sim.add_source('GenericSource', 'mybeam')
    source.particle = 'proton'
    source.energy.mono = energy * MeV
    source.energy.type = 'mono'
    source.position.type = 'point'
    source.position.translation = [0. * mm, 0. * mm, -.6 * mm]
    source.direction.type = 'momentum'
    source.direction.momentum = [0, 0, 1]

    source.n = 10 if sim.visu else number_of_particles

    # Physics list
    sim.physics_manager.physics_list_name = 'G4EmStandardPhysics_option4'

    # Phase spaces

    def add_detector(name, translation, attach_to_phantom=False):
        plane = sim.add_volume('Box', 'PlanePhaseSpace' + name)
        plane.size = [1. * mm, 1. * mm, 1 * nm]
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
            'Position',
            'KineticEnergy'
        ]
        particle_filter = sim.add_filter('ParticleFilter', 'Filter' + name)
        particle_filter.particle = 'proton'

        phase_space.filters.append(particle_filter)

    add_detector('In', [0. * mm, 0. * mm, -.5 * mm])
    add_detector('Out', [0. * mm, 0. * mm, .5 * mm])

    # Particle stats
    stat = sim.add_actor('SimulationStatisticsActor', 'stat')
    stat.output_filename = os.path.join(output_base, 'stats.txt')

    sim.run()


def spwatercalibration_fit(
    energy_max=DEFAULT_ENERGY_MAX,
    energy_min=DEFAULT_ENERGY_MIN,
    energy_samples=DEFAULT_ENERGY_SAMPLES,
    output=DEFAULT_OUTPUT,
    verbose=False
):
    def print_verbose(*args, **kwargs):
        if verbose:
            print(*args, **kwargs)

    sps = {}

    for energy in np.linspace(energy_max, energy_min, energy_samples):

        print_verbose("Processing events for energy of", energy, "MeV")
        output_base = output_format(output, energy)
        sps[energy] = []

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

                eloss = data['KineticEnergy'][event_mask][0] - data['KineticEnergy'][event_mask][1]

                sps[energy].append(eloss)

        except FileNotFoundError:
            continue

    def points_gen():
        for energy, elosses in sps.items():
            yield energy, np.median(elosses)

    points = np.array(list(points_gen()))
    x_e = points[:, 0]
    y_sp = points[:, 1]

    p = np.polyfit(x_e, y_sp, deg=10)
    print_verbose("Fitted parameters:", p)

    with open(os.path.join(output, 'fit.dat'), 'w', encoding='utf-8') as f:
        json.dump(list(p), f)
    with open(os.path.join(output, 'points.dat'), 'w', encoding='utf-8') as f:
        for x, y in zip(x_e, y_sp):
            f.write(f'{x} {y}\n')

    plt.figure()

    xs = np.linspace(energy_max, 0., 1000)
    plt.plot(xs, np.polyval(p, xs))

    plt.scatter(x_e, y_sp, marker='+')
    plt.xlabel("Energy [MeV]")
    plt.ylabel("Stopping power [MeV/mm]")

    plt.savefig(os.path.join(output, 'plot.pdf'))
    # plt.show()


def pctspwatercalibration(
    output=DEFAULT_OUTPUT,
    energy_max=DEFAULT_ENERGY_MAX,
    energy_min=DEFAULT_ENERGY_MIN,
    energy_samples=DEFAULT_ENERGY_SAMPLES,
    number_of_particles=DEFAULT_NUMBER_OF_PARTICLES,
    visu=False,
    verbose=False
):

    with Pool(maxtasksperchild=1) as pool:

        results = []

        energies = np.linspace(energy_max, energy_min, energy_samples)

        for energy in energies:
            result = pool.apply_async(pctwatercalibration_mc, kwds={
                'energy': energy,
                'number_of_particles': number_of_particles,
                'output': output,
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

    spwatercalibration_fit(
        energy_max,
        energy_min,
        energy_samples,
        output,
        verbose
    )


def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--output', help="Path of outputs", default=DEFAULT_OUTPUT)
    parser.add_argument('-e', '--energy-max', help="Maximum energy", default=DEFAULT_ENERGY_MAX, type=float)
    parser.add_argument('--energy-min', help="Minimum energy", default=DEFAULT_ENERGY_MIN, type=float)
    parser.add_argument('--energy-samples', help="Number of energy samples", default=DEFAULT_ENERGY_SAMPLES, type=int)
    parser.add_argument('-n', '--number-of-particles', help="Number of particles", default=DEFAULT_NUMBER_OF_PARTICLES, type=int)
    parser.add_argument('--visu', help="Visualize Monte Carlo simulation", default=False, action='store_true')
    parser.add_argument('--verbose', '-v', help="Verbose execution", default=False, action='store_true')
    args_info = parser.parse_args()

    pctspwatercalibration(**vars(args_info))

if __name__ == '__main__':
    main()
