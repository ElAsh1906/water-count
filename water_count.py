#!/usr/bin/python3.10

import pandas as pd
import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
import seaborn as sns
from MDAnalysis.analysis import distances
from tqdm import tqdm
import argparse


def count():
    """
    Count the number of water molecules between upper and lower halves of the membrane in each frame of the trajectory.

    Returns: resids (list): List containing residue IDs of water molecules between the upper and lower halves of the
    membrane in each frame.
    """
    resids = []
    with tqdm(total=num_frames, desc="Processing Frames") as pbar:
        for ts in u.trajectory:
            distance = distances.distance_array(upper_half, lower_half).min()
            water_between = mda.analysis.distances.between(molecule, upper_half, lower_half, distance - args.d)
            resids.append(water_between.resids)
            pbar.update(1)
    return resids


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculation of order parameters')
    parser.add_argument('-c', type=str, help='Topology file (default: %(default)s)', default='conf.gro',
                        metavar='<.gro>')
    parser.add_argument('-s', type=str, help='Topology file (default: %(default)s)', default='topol.tpr',
                        metavar='<.tpr>')
    parser.add_argument('-f', type=str, help='Trajectory file (default: %(default)s)', default='traj_comp.xtc',
                        metavar='<.xtc>')
    parser.add_argument('-o', type=str, help='Output text file (default: %(default)s)',
                        default='water_mols.csv', metavar='<.dat/.txt/.csv/...>')
    parser.add_argument('-l', '--lipids', type=str, help='Lipid model used, (default: %(default)s)',
                        default='DOPC', metavar='<lipid>')
    parser.add_argument('-la', type=str, nargs='+', help='Lipid atom used', metavar='<atom name>')
    parser.add_argument('-w', type=str, help='Resname of molecule to be counted, (default: %(default)s)',
                        default='SOL', metavar='<water>')
    parser.add_argument('-wa', type=str, help='Main atom of molecule to be counted, (default: %(default)s)',
                        default='OW', metavar='<atom>')
    parser.add_argument('-d', type=float, metavar='<distance>', default=0,
                        help='Distance for decreasing layer in which molecules are counted (default: %(default)sA)')
    parser.add_argument('-rw', type=int, metavar='<rolling window>', default=10,
                        help='Rolling window parameter for smoothed curve (default: %(default)s)')
    parser.add_argument('--save_resids', help='Save all resids', action='store_true')
    parser.add_argument('--calc', help='Calculate mean and error', action='store_true')
    parser.add_argument('--plot', help='Plot data', action='store_true')

    args = parser.parse_args()

    if args.s != 'topol.tpr':
        u = mda.Universe(args.s, args.f)
    else:
        u = mda.Universe(args.c, args.f)

    num_frames = len(u.trajectory)
    time = u.trajectory.totaltime

    membrane = u.select_atoms(f"resname {args.lipids}")
    membrane_com = membrane.center_of_mass(compound='group')
    z_cutoff = membrane_com[2]

    la = ' or '.join(['name ' + i for i in args.la])
    lipid_atoms = f"({la} and resname {args.lipids})"

    upper_half = u.select_atoms(f"{lipid_atoms} and (prop z > {z_cutoff})")
    lower_half = u.select_atoms(f"{lipid_atoms} and ( prop z < {z_cutoff})")

    molecule = u.select_atoms(f"resname {args.w} and name {args.wa}")

    final_arr = count()
    sum_arr = np.array([len(subarray) for subarray in final_arr])

    np.savetxt(args.o, sum_arr, delimiter=',', fmt='%.4f')

    if args.save_resids:
        with open("resids.csv", "w") as f:
            for row in final_arr:
                f.write("%s\n" % ','.join(str(col) for col in row))

    if args.calc:
        print(sum_arr.mean().round(2), u"\u00B1", round(sum_arr.std() / np.sqrt(num_frames), 2))

    if args.plot:
        smoothed = pd.DataFrame(sum_arr).rolling(args.rw).mean()
        plt.figure(figsize=(12, 5))
        plt.xlabel('t, ns')

        ax1 = sns.lineplot(y=sum_arr, x=np.linspace(0, time, len(sum_arr)) / 1000, color='blue',
                           label='Number of water')
        ax2 = sns.lineplot(y=smoothed[0], x=np.linspace(0, time, len(smoothed[0])) / 1000, color='red', label='Rolling '
                                                                                                              'window')

        plt.legend(frameon=False)
        plt.show()
