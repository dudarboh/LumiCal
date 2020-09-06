from ROOT import TFile, TTree

import time
import numpy as np
from clustering import make_clusters_lists

from output_tree import OutputTree
from hit import Hit, TrHit
import argparse


def extract_noise():
    noise = np.zeros((4, 64, 8)) # In MeV
    f=open("./noise.txt", "r")
    lines = f.readlines()
    for line in lines[1:]:
        values = line.split('  ')
        sector = int(values[0])
        if sector == -1:
            continue
        pad = int(values[1])
        layer = int(values[2])
        noise[sector, pad, layer] = float(values[3])
    return noise


def bad_pad(sector, pad, layer):
    """Return true if pad is bad"""
    return ((layer == 0 and sector == 1 and pad in (62,))
            or (layer == 0 and sector == 2 and pad in (20, 22, 57, 61, 63))
            or (layer == 2 and sector == 1 and pad in (28, 31, 34, 63))
            or (layer == 2 and sector == 2 and pad in (38, 53, 62))
            or (layer == 3 and sector == 1 and pad in (63,))
            or (layer == 3 and sector == 2 and pad in (31, 33, 52, 55, 61, 62))
            or (layer == 4 and sector == 1 and pad in (29, 39, 41, 55, 56, 63))
            or (layer == 4 and sector == 2 and pad in (28, 62))
            or (layer == 5 and sector == 1 and pad in (32, 36, 40, 41, 44, 45, 49, 56, 58, 63))
            or (layer == 5 and sector == 2 and pad in (28, 52, 54, 61, 62))
            or (layer == 6 and sector == 1 and pad in (26, 30, 62, 63))
            or (layer == 6 and sector == 2 and pad in (34, 42, 54, 57, 59, 60, 62))
            or sector == -1)  # grounded channel




def make_hits_lists(event):
    hits_tracker1 = []
    hits_tracker2 = []
    hits_calorimeter = []



    for (sector, pad, layer, energy,
         hit_type, track_len,
         p_x, p_y, p_z, p_px, p_py, p_pz, p_e) in zip(event.tr1_sector,
                                                      event.tr1_pad,
                                                      event.tr1_layer,
                                                      event.tr1_energy,
                                                      event.tr1_type,
                                                      event.tr1_track_len,
                                                      event.tr1_x,
                                                      event.tr1_y,
                                                      event.tr1_z,
                                                      event.tr1_px,
                                                      event.tr1_py,
                                                      event.tr1_pz,
                                                      event.tr1_p_energy):

        # Selection as in data for tracker1
        if (pad < 20 or sector == 0 or sector == 3
            or energy <= 0. or bad_pad(sector, pad, layer)):
            continue

        # If passed the selection create hit and add to corresponding list
        hit = TrHit(sector, pad, layer, energy,
               hit_type, track_len,
               p_x, p_y, p_z, p_px, p_py, p_pz, p_e)

        hits_tracker1.append(hit)

    for (sector, pad, layer, energy,
         hit_type, track_len,
         p_x, p_y, p_z, p_px, p_py, p_pz, p_e) in zip(event.tr2_sector,
                                                      event.tr2_pad,
                                                      event.tr2_layer,
                                                      event.tr2_energy,
                                                      event.tr2_type,
                                                      event.tr2_track_len,
                                                      event.tr2_x,
                                                      event.tr2_y,
                                                      event.tr2_z,
                                                      event.tr2_px,
                                                      event.tr2_py,
                                                      event.tr2_pz,
                                                      event.tr2_p_energy):
        # Selection as in data for tracker2
        if (pad < 20 or sector == 0 or sector == 3
            or energy <= 0. or bad_pad(sector, pad, layer)):
            continue

        # If passed the selection create hit and add to corresponding list
        hits_tracker2.append(TrHit(sector, pad, layer, energy,
                                   hit_type, track_len,
                                   p_x, p_y, p_z, p_px, p_py, p_pz, p_e))

    for (sector, pad, layer, energy) in zip(event.cal_sector,
                                            event.cal_pad,
                                            event.cal_layer,
                                            event.cal_energy):

        # Selection as in data for calorimeter
        if (pad < 20 or layer == 7
            or sector == 0 or sector == 3
            or bad_pad(sector, pad, layer)):
            continue

        # If passed the selection create hit and add to corresponding list
        hits_calorimeter.append(Hit(sector, pad, layer, energy))

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
    input_tree = input_file.lumical
    print("Total n events in loaded files: ", input_tree.GetEntries())

    # Create output root file.
    # Create output root file before the tree!!! It prevents memory leakage
    output_file = OutputTree()

    noise = extract_noise()

    n_events = input_tree.GetEntries()
    for idx, event in enumerate(input_tree):
        if idx == 20000:
            break

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

        output_file.fill_output_tree(event, hits_tr1, hits_tr2, hits_cal, clusters_tr1, clusters_tr2, clusters_cal)

    output_file.write_file()
    print("Hooray, extracted tree file is ready, take it :3")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=('Do selection and clustering of data/MC'))
    parser.add_argument('path_to_file', type=str, help='Provide path to the root file')
    args = parser.parse_args()

    main(args)
