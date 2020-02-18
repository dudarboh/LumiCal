from ROOT import TFile, TTree
import array
import time
import numpy as np
from clustering import make_clusters_lists

import cProfile
import pstats

import argparse


def bad_pad(sector, pad, layer):
    return ((layer == 0 and sector == 1 and pad in (26, 61))
            or (layer == 0 and sector == 2 and pad in (31, 57, 61))
            or (layer == 2 and sector == 1 and pad in (28, 31, 34))
            or (layer == 2 and sector == 2 and pad in (38, 53))
            or (layer == 3 and sector == 2 and pad in (31, 33, 52, 55, 61))
            or (layer == 4 and sector == 1 and pad in (29, 39, 41, 55, 56))
            or (layer == 4 and sector == 2 and pad in (28,))
            or (layer == 5 and sector == 1 and pad in (32, 36, 40, 41, 44, 45, 49, 56, 58))
            or (layer == 5 and sector == 2 and pad in (28, 52, 54, 61))
            or (layer == 6 and sector == 1 and pad in (26, 30))
            or (layer == 6 and sector == 2 and pad in (34, 42, 54, 57, 59, 60)))


class TrHit:
    def __init__(self, sector, pad, layer, energy_in_mip,
                 hit_type, track_len, track_id,
                 p_x, p_y, p_z, p_px, p_py, p_pz, p_energy):
        self.sector = sector
        self.pad = pad
        self.layer = layer
        self.energy = energy_in_mip

        self.type = hit_type
        self.track_len = track_len
        self.track_id = track_id

        self.p_x = p_x
        self.p_y = p_y
        self.p_z = p_z
        self.p_px = p_px
        self.p_py = p_py
        self.p_pz = p_pz
        self.p_energy = p_energy
        rho = 80. + 0.9 + 1.8 * self.pad
        phi = np.pi / 2 + np.pi / 12 - np.pi / 48 - np.pi / 24 * self.sector
        self.x = rho * np.cos(phi)
        self.y = rho * np.sin(phi)
        self.seed = -1


class CalHit:
    def __init__(self, sector, pad, layer, energy_in_mip):
        self.sector = sector
        self.pad = pad
        self.layer = layer
        self.energy = energy_in_mip
        rho = 80. + 0.9 + 1.8 * self.pad
        phi = np.pi / 2 + np.pi / 12 - np.pi / 48 - np.pi / 24 * self.sector
        self.x = rho * np.cos(phi)
        self.y = rho * np.sin(phi)


def make_hits_lists(event):
    hits_tracker1 = []
    hits_tracker2 = []
    hits_calorimeter = []

    for (sector, pad, layer, energy,
         hit_type, track_len, track_id,
         p_x, p_y, p_z, p_px, p_py, p_pz, p_e) in zip(event.tr1_hit_sector,
                                                      event.tr1_hit_pad,
                                                      event.tr1_hit_layer,
                                                      event.tr1_hit_energy,
                                                      event.tr1_hit_type,
                                                      event.tr1_hit_track_len,
                                                      event.tr1_hit_track_id,
                                                      event.tr1_particle_x,
                                                      event.tr1_particle_y,
                                                      event.tr1_particle_z,
                                                      event.tr1_particle_px,
                                                      event.tr1_particle_py,
                                                      event.tr1_particle_pz,
                                                      event.tr1_particle_energy):

        # Selection as in tracker1
        if (pad < 20 or sector == 0 or sector == 3
            or energy <= 0. or bad_pad(sector, pad, layer)):
            continue

        # If passed the selection create hit and add to corresponding list
        hits_tracker1.append(TrHit(sector, pad, layer, energy,
                                   hit_type, track_len, track_id,
                                   p_x, p_y, p_z, p_px, p_py, p_pz, p_e))

    for (sector, pad, layer, energy,
         hit_type, track_len, track_id,
         p_x, p_y, p_z, p_px, p_py, p_pz, p_e) in zip(event.tr2_hit_sector,
                                                      event.tr2_hit_pad,
                                                      event.tr2_hit_layer,
                                                      event.tr2_hit_energy,
                                                      event.tr2_hit_type,
                                                      event.tr2_hit_track_len,
                                                      event.tr2_hit_track_id,
                                                      event.tr2_particle_x,
                                                      event.tr2_particle_y,
                                                      event.tr2_particle_z,
                                                      event.tr2_particle_px,
                                                      event.tr2_particle_py,
                                                      event.tr2_particle_pz,
                                                      event.tr2_particle_energy):
        # Selection as in data
        if (pad < 20 or sector == 0 or sector == 3
            or energy <= 0. or bad_pad(sector, pad, layer)):
            continue

        # If passed the selection create hit and add to corresponding list
        hits_tracker2.append(TrHit(sector, pad, layer, energy,
                                   hit_type, track_len, track_id,
                                   p_x, p_y, p_z, p_px, p_py, p_pz, p_e))

    for (sector, pad, layer, energy) in zip(event.cal_hit_sector,
                                            event.cal_hit_pad,
                                            event.cal_hit_layer,
                                            event.cal_hit_energy):

        # Selection as in data for calorimeter
        if (energy < 1.4
            or pad < 20 or layer == 7
            or sector == 0 or sector == 3
            or bad_pad(sector, pad, layer)):
            continue

        # If passed the selection create hit and add to corresponding list
        hits_calorimeter.append(CalHit(sector, pad, layer, energy))

    return hits_tracker1, hits_tracker2, hits_calorimeter


def align_mc(hits_tr1, hits_tr2, hits_cal):
    """Shift hit's y coordinates for misalignment"""
    tr1_shift = -0.2480096316087952
    tr2_shift = 0.9012092661162399
    cal_shift = -0.6531996345075299

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
    input_tree = input_file.mc
    print("Total n events in loaded files: ", input_tree.GetEntries())

    # Create output root file.
    # Create output root file before the tree!!! It prevents memory leakage
    output_file = TFile("./extracted_mc_RENAME.root", "RECREATE")
    output_tree = TTree('mc', 'Extracted MC')

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

    tr1_hit_type = array.array('i', [0] * 128)
    tr1_track_len = array.array('f', [0.0] * 128)
    tr1_track_id = array.array('i', [0] * 128)
    tr1_particle_x = array.array('f', [0.0] * 128)
    tr1_particle_y = array.array('f', [0.0] * 128)
    tr1_particle_z = array.array('f', [0.0] * 128)
    tr1_particle_px = array.array('f', [0.0] * 128)
    tr1_particle_py = array.array('f', [0.0] * 128)
    tr1_particle_pz = array.array('f', [0.0] * 128)
    tr1_particle_energy = array.array('f', [0.0] * 128)

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

    tr2_hit_type = array.array('i', [0] * 128)
    tr2_track_len = array.array('f', [0.0] * 128)
    tr2_track_id = array.array('i', [0] * 128)
    tr2_particle_x = array.array('f', [0.0] * 128)
    tr2_particle_y = array.array('f', [0.0] * 128)
    tr2_particle_z = array.array('f', [0.0] * 128)
    tr2_particle_px = array.array('f', [0.0] * 128)
    tr2_particle_py = array.array('f', [0.0] * 128)
    tr2_particle_pz = array.array('f', [0.0] * 128)
    tr2_particle_energy = array.array('f', [0.0] * 128)

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

    output_tree.Branch('tr1_hit_type', tr1_hit_type, 'tr1_hit_type[tr1_n_hits]/I')
    output_tree.Branch('tr1_track_len', tr1_track_len, 'tr1_track_len[tr1_n_hits]/F')
    output_tree.Branch('tr1_track_id', tr1_track_id, 'tr1_track_id[tr1_n_hits]/I')
    output_tree.Branch('tr1_particle_x', tr1_particle_x, 'tr1_particle_x[tr1_n_hits]/F')
    output_tree.Branch('tr1_particle_y', tr1_particle_y, 'tr1_particle_y[tr1_n_hits]/F')
    output_tree.Branch('tr1_particle_z', tr1_particle_z, 'tr1_particle_z[tr1_n_hits]/F')
    output_tree.Branch('tr1_particle_px', tr1_particle_px, 'tr1_particle_px[tr1_n_hits]/F')
    output_tree.Branch('tr1_particle_py', tr1_particle_py, 'tr1_particle_py[tr1_n_hits]/F')
    output_tree.Branch('tr1_particle_pz', tr1_particle_pz, 'tr1_particle_pz[tr1_n_hits]/F')
    output_tree.Branch('tr1_particle_energy', tr1_particle_energy, 'tr1_particle_energy[tr1_n_hits]/F')

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

    output_tree.Branch('tr2_hit_type', tr2_hit_type, 'tr2_hit_type[tr2_n_hits]/I')
    output_tree.Branch('tr2_track_len', tr2_track_len, 'tr2_track_len[tr2_n_hits]/F')
    output_tree.Branch('tr2_track_id', tr2_track_id, 'tr2_track_id[tr2_n_hits]/I')
    output_tree.Branch('tr2_particle_x', tr2_particle_x, 'tr2_particle_x[tr2_n_hits]/F')
    output_tree.Branch('tr2_particle_y', tr2_particle_y, 'tr2_particle_y[tr2_n_hits]/F')
    output_tree.Branch('tr2_particle_z', tr2_particle_z, 'tr2_particle_z[tr2_n_hits]/F')
    output_tree.Branch('tr2_particle_px', tr2_particle_px, 'tr2_particle_px[tr2_n_hits]/F')
    output_tree.Branch('tr2_particle_py', tr2_particle_py, 'tr2_particle_py[tr2_n_hits]/F')
    output_tree.Branch('tr2_particle_pz', tr2_particle_pz, 'tr2_particle_pz[tr2_n_hits]/F')
    output_tree.Branch('tr2_particle_energy', tr2_particle_energy, 'tr2_particle_energy[tr2_n_hits]/F')

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
        # if idx == 1000:
        #     break

        if idx % (1000) == 0:
            time_min = (time.time() - start_time) // 60
            time_sec = (time.time() - start_time) % 60
            print('Event: {} out of {};'.format(idx, n_events), end=' ')
            print('{} min {} sec'.format(time_min, time_sec))

        hits_tr1, hits_tr2, hits_cal = make_hits_lists(event)
        align_mc(hits_tr1, hits_tr2, hits_cal)
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

            tr1_hit_type[i] = hit.type
            tr1_track_len[i] = hit.track_len
            tr1_track_id[i] = hit.track_id
            tr1_particle_x[i] = hit.p_x
            tr1_particle_y[i] = hit.p_y
            tr1_particle_z[i] = hit.p_z
            tr1_particle_px[i] = hit.p_px
            tr1_particle_py[i] = hit.p_py
            tr1_particle_pz[i] = hit.p_pz
            tr1_particle_energy[i] = hit.p_energy

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

            tr2_hit_type[i] = hit.type
            tr2_track_len[i] = hit.track_len
            tr2_track_id[i] = hit.track_id
            tr2_particle_x[i] = hit.p_x
            tr2_particle_y[i] = hit.p_y
            tr2_particle_z[i] = hit.p_z
            tr2_particle_px[i] = hit.p_px
            tr2_particle_py[i] = hit.p_py
            tr2_particle_pz[i] = hit.p_pz
            tr2_particle_energy[i] = hit.p_energy

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

    pr = cProfile.Profile()
    pr.enable()
    main(args)
    pr.disable()

    ps = pstats.Stats(pr).sort_stats(pstats.SortKey.CUMULATIVE)
    ps.print_stats()
