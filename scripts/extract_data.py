from ROOT import TFile, TTree
import array
import time
import numpy as np
from clustering import make_clusters_lists

import argparse


class Hit:
    def __init__(self, sector, pad, layer, energy):
        self.sector = sector
        self.pad = pad
        self.layer = layer
        self.energy = energy

        rho = 80.9 + 1.8 * pad
        phi = np.pi / 24. * (13.5 - sector)
        self.x = rho * np.cos(phi)
        self.y = rho * np.sin(phi)

        if layer == 0 or layer == 1:
            self.seed = -1


def make_hits_lists(event):
    """Return hits lists from a root file event"""

    hits_calorimeter = []
    hits_tracker1 = []
    hits_tracker2 = []

    for sector, pad, layer, energy in zip(event.sector, event.pad,
                                          event.layer, event.energy):

        # Geometrical cuts
        if (pad < 20 or sector == 0 or sector == 3 or layer == 7
           or (layer > 1 and energy <= 1.4)):
            continue

        if layer == 0:
            hits_tracker1.append(Hit(sector, pad, layer, energy))
        elif layer == 1:
            hits_tracker2.append(Hit(sector, pad, layer, energy))
        else:
            hits_calorimeter.append(Hit(sector, pad, layer, energy))

    return hits_tracker1, hits_tracker2, hits_calorimeter


def align_data(hits_tr1, hits_tr2, hits_cal):
    """Shift hit's y coordinates for misalignment"""
    tr1_shift = -0.14156476251841354
    tr2_shift = 0.9273328597379873
    cal_shift = -0.785768097219659

    for hit in hits_tr1:
        hit.y -= tr1_shift
    for hit in hits_tr2:
        hit.y -= tr2_shift
    for hit in hits_cal:
        hit.y -= cal_shift


def main(args):
    start_time = time.time()

    # Upload data for analysis
    input_file = TFile.Open(args.path_to_file, "READ")
    input_tree = input_file.data
    print("Total n events in loaded files: ", input_tree.GetEntries())

    # Create output root file.
    # Create output root file before the tree!!! It prevents memory leakage
    output_file = TFile('./extracted_data_RENAME.root', "RECREATE")
    output_tree = TTree('data', 'Extracted Data')

    # Create variables associated with the output tree
    tr1_n_hits = array.array('i', [0])
    tr1_hit_pad = array.array('i', [0] * 128)
    tr1_hit_sector = array.array('i', [0] * 128)
    tr1_hit_layer = array.array('i', [0] * 128)
    tr1_hit_x = array.array('f', [0.0] * 128)
    tr1_hit_y = array.array('f', [0.0] * 128)
    tr1_hit_energy = array.array('f', [0.0] * 128)
    tr1_n_clusters = array.array('i', [0])
    tr1_cluster_n_pads = array.array('i', [0] * 128)
    tr1_cluster_pad = array.array('f', [0.0] * 128)
    tr1_cluster_sector = array.array('f', [0.0] * 128)
    tr1_cluster_x = array.array('f', [0.0] * 128)
    tr1_cluster_y = array.array('f', [0.0] * 128)
    tr1_cluster_energy = array.array('f', [0.0] * 128)

    tr2_n_hits = array.array('i', [0])
    tr2_hit_pad = array.array('i', [0] * 128)
    tr2_hit_sector = array.array('i', [0] * 128)
    tr2_hit_layer = array.array('i', [0] * 128)
    tr2_hit_x = array.array('f', [0.0] * 128)
    tr2_hit_y = array.array('f', [0.0] * 128)
    tr2_hit_energy = array.array('f', [0.0] * 128)
    tr2_n_clusters = array.array('i', [0])
    tr2_cluster_n_pads = array.array('i', [0] * 128)
    tr2_cluster_pad = array.array('f', [0.0] * 128)
    tr2_cluster_sector = array.array('f', [0.0] * 128)
    tr2_cluster_x = array.array('f', [0.0] * 128)
    tr2_cluster_y = array.array('f', [0.0] * 128)
    tr2_cluster_energy = array.array('f', [0.0] * 128)

    cal_n_hits = array.array('i', [0])
    cal_hit_pad = array.array('i', [0] * 128 * 5)
    cal_hit_sector = array.array('i', [0] * 128 * 5)
    cal_hit_layer = array.array('i', [0] * 128 * 5)
    cal_hit_x = array.array('f', [0.0] * 128 * 5)
    cal_hit_y = array.array('f', [0.0] * 128 * 5)
    cal_hit_energy = array.array('f', [0.0] * 128 * 5)
    cal_n_clusters = array.array('i', [0])
    cal_cluster_n_pads = array.array('i', [0] * 128 * 5)
    cal_cluster_n_towers = array.array('i', [0] * 128 * 5)
    cal_cluster_pad = array.array('f', [0.0] * 128 * 5)
    cal_cluster_sector = array.array('f', [0.0] * 128 * 5)
    cal_cluster_layer = array.array('f', [0.0] * 128 * 5)
    cal_cluster_x = array.array('f', [0.0] * 128 * 5)
    cal_cluster_y = array.array('f', [0.0] * 128 * 5)
    cal_cluster_energy = array.array('f', [0.0] * 128 * 5)

    # Create branches in the output tree for these variables
    output_tree.Branch('tr1_n_hits', tr1_n_hits, 'tr1_n_hits/I')
    output_tree.Branch('tr1_hit_pad', tr1_hit_pad, 'tr1_hit_pad[tr1_n_hits]/I')
    output_tree.Branch('tr1_hit_sector', tr1_hit_sector, 'tr1_hit_sector[tr1_n_hits]/I')
    output_tree.Branch('tr1_hit_layer', tr1_hit_layer, 'tr1_hit_layer[tr1_n_hits]/I')
    output_tree.Branch('tr1_hit_x', tr1_hit_x, 'tr1_hit_x[tr1_n_hits]/F')
    output_tree.Branch('tr1_hit_y', tr1_hit_y, 'tr1_hit_y[tr1_n_hits]/F')
    output_tree.Branch('tr1_hit_energy', tr1_hit_energy, 'tr1_hit_energy[tr1_n_hits]/F')
    output_tree.Branch('tr1_n_clusters', tr1_n_clusters, 'tr1_n_clusters/I')
    output_tree.Branch('tr1_cluster_n_pads', tr1_cluster_n_pads, 'tr1_cluster_n_pads[tr1_n_clusters]/I')
    output_tree.Branch('tr1_cluster_pad', tr1_cluster_pad, 'tr1_cluster_pad[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_sector', tr1_cluster_sector, 'tr1_cluster_sector[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_x', tr1_cluster_x, 'tr1_cluster_x[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_y', tr1_cluster_y, 'tr1_cluster_y[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_energy', tr1_cluster_energy, 'tr1_cluster_energy[tr1_n_clusters]/F')

    output_tree.Branch('tr2_n_hits', tr2_n_hits, 'tr2_n_hits/I')
    output_tree.Branch('tr2_hit_pad', tr2_hit_pad, 'tr2_hit_pad[tr2_n_hits]/I')
    output_tree.Branch('tr2_hit_sector', tr2_hit_sector, 'tr2_hit_sector[tr2_n_hits]/I')
    output_tree.Branch('tr2_hit_layer', tr2_hit_layer, 'tr2_hit_layer[tr2_n_hits]/I')
    output_tree.Branch('tr2_hit_x', tr2_hit_x, 'tr2_hit_x[tr2_n_hits]/F')
    output_tree.Branch('tr2_hit_y', tr2_hit_y, 'tr2_hit_y[tr2_n_hits]/F')
    output_tree.Branch('tr2_hit_energy', tr2_hit_energy, 'tr2_hit_energy[tr2_n_hits]/F')
    output_tree.Branch('tr2_n_clusters', tr2_n_clusters, 'tr2_n_clusters/I')
    output_tree.Branch('tr2_cluster_n_pads', tr2_cluster_n_pads, 'tr2_cluster_n_pads[tr2_n_clusters]/I')
    output_tree.Branch('tr2_cluster_pad', tr2_cluster_pad, 'tr2_cluster_pad[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_sector', tr2_cluster_sector, 'tr2_cluster_sector[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_x', tr2_cluster_x, 'tr2_cluster_x[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_y', tr2_cluster_y, 'tr2_cluster_y[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_energy', tr2_cluster_energy, 'tr2_cluster_energy[tr2_n_clusters]/F')

    output_tree.Branch('cal_n_hits', cal_n_hits, 'cal_n_hits/I')
    output_tree.Branch('cal_hit_pad', cal_hit_pad, 'cal_hit_pad[cal_n_hits]/I')
    output_tree.Branch('cal_hit_sector', cal_hit_sector, 'cal_hit_sector[cal_n_hits]/I')
    output_tree.Branch('cal_hit_layer', cal_hit_layer, 'cal_hit_layer[cal_n_hits]/I')
    output_tree.Branch('cal_hit_x', cal_hit_x, 'cal_hit_x[cal_n_hits]/F')
    output_tree.Branch('cal_hit_y', cal_hit_y, 'cal_hit_y[cal_n_hits]/F')
    output_tree.Branch('cal_hit_energy', cal_hit_energy, 'cal_hit_energy[cal_n_hits]/F')
    output_tree.Branch('cal_n_clusters', cal_n_clusters, 'cal_n_clusters/I')
    output_tree.Branch('cal_cluster_n_pads', cal_cluster_n_pads, 'cal_cluster_n_pads[cal_n_clusters]/I')
    output_tree.Branch('cal_cluster_n_towers', cal_cluster_n_towers, 'cal_cluster_n_towers[cal_n_clusters]/I')
    output_tree.Branch('cal_cluster_pad', cal_cluster_pad, 'cal_cluster_pad[cal_n_clusters]/F')
    output_tree.Branch('cal_cluster_sector', cal_cluster_sector, 'cal_cluster_sector[cal_n_clusters]/F')
    output_tree.Branch('cal_cluster_layer', cal_cluster_layer, 'cal_cluster_layer[cal_n_clusters]/F')
    output_tree.Branch('cal_cluster_x', cal_cluster_x, 'cal_cluster_x[cal_n_clusters]/F')
    output_tree.Branch('cal_cluster_y', cal_cluster_y, 'cal_cluster_y[cal_n_clusters]/F')
    output_tree.Branch('cal_cluster_energy', cal_cluster_energy, 'cal_cluster_energy[cal_n_clusters]/F')

    n_events = input_tree.GetEntries()
    for idx, event in enumerate(input_tree):
        if idx != 5:
            continue

        if idx % (1000) == 0:
            time_min = (time.time() - start_time) // 60
            time_sec = (time.time() - start_time) % 60
            print('Event: {} out of {};'.format(idx, n_events), end=' ')
            print('{} min {} sec'.format(time_min, time_sec))

        hits_tr1, hits_tr2, hits_cal = make_hits_lists(event)
        align_data(hits_tr1, hits_tr2, hits_cal)
        clusters_tr1, clusters_tr2, clusters_cal = make_clusters_lists(hits_tr1, hits_tr2, hits_cal)

        # Resort clusters in trackers by distance to main cluster in calorimeter
        if len(clusters_cal) != 0:
            clusters_tr1.sort(key=lambda x: abs(x.y - clusters_cal[0].y))
            clusters_tr2.sort(key=lambda x: abs(x.y - clusters_cal[0].y))

        tr1_n_hits[0] = len(hits_tr1)
        for i, hit in enumerate(hits_tr1):
            tr1_hit_pad[i] = hit.pad
            tr1_hit_sector[i] = hit.sector
            tr1_hit_layer[i] = hit.layer
            tr1_hit_energy[i] = hit.energy
            tr1_hit_x[i] = hit.x
            tr1_hit_y[i] = hit.y
        tr1_n_clusters[0] = len(clusters_tr1)
        for i, cluster in enumerate(clusters_tr1):
            tr1_cluster_n_pads[i] = cluster.n_pads
            tr1_cluster_pad[i] = cluster.pad
            tr1_cluster_sector[i] = cluster.sector
            tr1_cluster_x[i] = cluster.x
            tr1_cluster_y[i] = cluster.y
            tr1_cluster_energy[i] = cluster.energy

        tr2_n_hits[0] = len(hits_tr2)
        for i, hit in enumerate(hits_tr2):
            tr2_hit_pad[i] = hit.pad
            tr2_hit_sector[i] = hit.sector
            tr2_hit_layer[i] = hit.layer
            tr2_hit_energy[i] = hit.energy
            tr2_hit_x[i] = hit.x
            tr2_hit_y[i] = hit.y
        tr2_n_clusters[0] = len(clusters_tr2)
        for i, cluster in enumerate(clusters_tr2):
            tr2_cluster_n_pads[i] = cluster.n_pads
            tr2_cluster_pad[i] = cluster.pad
            tr2_cluster_sector[i] = cluster.sector
            tr2_cluster_x[i] = cluster.x
            tr2_cluster_y[i] = cluster.y
            tr2_cluster_energy[i] = cluster.energy

        cal_n_hits[0] = len(hits_cal)
        for i, hit in enumerate(hits_cal):
            cal_hit_pad[i] = hit.pad
            cal_hit_sector[i] = hit.sector
            cal_hit_layer[i] = hit.layer
            cal_hit_x[i] = hit.x
            cal_hit_y[i] = hit.y
            cal_hit_energy[i] = hit.energy
        cal_n_clusters[0] = len(clusters_cal)
        for i, cluster in enumerate(clusters_cal):
            cal_cluster_n_pads[i] = cluster.n_pads
            cal_cluster_n_towers[i] = cluster.n_towers
            cal_cluster_pad[i] = cluster.pad
            cal_cluster_sector[i] = cluster.sector
            cal_cluster_layer[i] = cluster.layer
            cal_cluster_x[i] = cluster.x
            cal_cluster_y[i] = cluster.y
            cal_cluster_energy[i] = cluster.energy

        output_tree.Fill()

    output_tree.Write()
    output_file.Close()

    print("Hooray, extracted tree file is ready, take it :3")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=('Do selection and clustering of data/MC'))
    parser.add_argument('path_to_file', type=str, help='Provide path to the root file')
    args = parser.parse_args()

    main(args)