from ROOT import TFile, TTree
import array
import time
import numpy as np

import argparse

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

class Hit:
    def __init__(self, sector, pad, layer, energy_in_mev):
        self.sector = sector
        self.pad = pad
        self.layer = layer
        self.energy = energy_in_mev

        rho = 80. + 0.9 + 1.8 * self.pad
        phi = np.pi / 2 + np.pi / 12 - np.pi / 48 - np.pi / 24 * self.sector
        self.x = rho * np.cos(phi)
        self.y = rho * np.sin(phi)

def make_hits_lists(event):
    hits_tracker1 = []
    hits_tracker2 = []
    hits_calorimeter = []

    for (sector, pad, layer, energy) in zip(event.tr1_sector,
                                            event.tr1_pad,
                                            event.tr1_layer,
                                            event.tr1_energy):

        # Selection as in data for tracker1
        if (pad < 20 or sector == 0 or sector == 3
            or energy <= 0. or bad_pad(sector, pad, layer)):
            continue

        # If passed the selection create hit and add to corresponding list
        hits_tracker1.append(Hit(sector, pad, layer, energy))

    for (sector, pad, layer, energy) in zip(event.tr2_sector,
                                           event.tr2_pad,
                                           event.tr2_layer,
                                           event.tr2_energy):
        # Selection as in data for tracker2
        if (pad < 20 or sector == 0 or sector == 3
            or energy <= 0. or bad_pad(sector, pad, layer)):
            continue

        # If passed the selection create hit and add to corresponding list
        hits_tracker2.append(Hit(sector, pad, layer, energy))

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

def main(args):
    start_time = time.time()
    # Upload data for analysis
    input_file = TFile.Open(args.path_to_file, "READ")
    input_tree = input_file.lumical
    print("Total n events in loaded files: ", input_tree.GetEntries())

    # Create output root file.
    # Create output root file before the tree!!! It prevents memory leakage
    output_file = TFile("./lucas_geocuts.root", "RECREATE")
    output_tree = TTree('lumical', 'MC')

    # Create variables associated with the output tree
    n_triggers = array.array('i', [0])
    trigger1 = array.array('i', [0])
    trigger2 = array.array('i', [0])
    trigger3 = array.array('i', [0])

    tr1_n_hits = array.array('i', [0])
    tr1_hit_pad = array.array('i', [0] * 128)
    tr1_hit_sector = array.array('i', [0] * 128)
    tr1_hit_layer = array.array('i', [0] * 128)
    tr1_hit_x = array.array('f', [0.0] * 128)
    tr1_hit_y = array.array('f', [0.0] * 128)
    tr1_hit_energy = array.array('f', [0.0] * 128)

    tr2_n_hits = array.array('i', [0])
    tr2_hit_pad = array.array('i', [0] * 128)
    tr2_hit_sector = array.array('i', [0] * 128)
    tr2_hit_layer = array.array('i', [0] * 128)
    tr2_hit_x = array.array('f', [0.0] * 128)
    tr2_hit_y = array.array('f', [0.0] * 128)
    tr2_hit_energy = array.array('f', [0.0] * 128)

    cal_n_hits = array.array('i', [0])
    cal_hit_pad = array.array('i', [0] * 128 * 5)
    cal_hit_sector = array.array('i', [0] * 128 * 5)
    cal_hit_layer = array.array('i', [0] * 128 * 5)
    cal_hit_x = array.array('f', [0.0] * 128 * 5)
    cal_hit_y = array.array('f', [0.0] * 128 * 5)
    cal_hit_energy = array.array('f', [0.0] * 128 * 5)

    # Create branches in the output tree for these variables
    output_tree.Branch('n_triggers', n_triggers, 'n_triggers/I')
    output_tree.Branch('trigger1', trigger1, 'trigger1/I')
    output_tree.Branch('trigger2', trigger2, 'trigger2/I')
    output_tree.Branch('trigger3', trigger3, 'trigger3/I')

    output_tree.Branch('tr1_n_hits', tr1_n_hits, 'tr1_n_hits/I')
    output_tree.Branch('tr1_hit_pad', tr1_hit_pad, 'tr1_hit_pad[tr1_n_hits]/I')
    output_tree.Branch('tr1_hit_sector', tr1_hit_sector, 'tr1_hit_sector[tr1_n_hits]/I')
    output_tree.Branch('tr1_hit_layer', tr1_hit_layer, 'tr1_hit_layer[tr1_n_hits]/I')
    output_tree.Branch('tr1_hit_x', tr1_hit_x, 'tr1_hit_x[tr1_n_hits]/F')
    output_tree.Branch('tr1_hit_y', tr1_hit_y, 'tr1_hit_y[tr1_n_hits]/F')
    output_tree.Branch('tr1_hit_energy', tr1_hit_energy, 'tr1_hit_energy[tr1_n_hits]/F')

    output_tree.Branch('tr2_n_hits', tr2_n_hits, 'tr2_n_hits/I')
    output_tree.Branch('tr2_hit_pad', tr2_hit_pad, 'tr2_hit_pad[tr2_n_hits]/I')
    output_tree.Branch('tr2_hit_sector', tr2_hit_sector, 'tr2_hit_sector[tr2_n_hits]/I')
    output_tree.Branch('tr2_hit_layer', tr2_hit_layer, 'tr2_hit_layer[tr2_n_hits]/I')
    output_tree.Branch('tr2_hit_x', tr2_hit_x, 'tr2_hit_x[tr2_n_hits]/F')
    output_tree.Branch('tr2_hit_y', tr2_hit_y, 'tr2_hit_y[tr2_n_hits]/F')
    output_tree.Branch('tr2_hit_energy', tr2_hit_energy, 'tr2_hit_energy[tr2_n_hits]/F')

    output_tree.Branch('cal_n_hits', cal_n_hits, 'cal_n_hits/I')
    output_tree.Branch('cal_hit_pad', cal_hit_pad, 'cal_hit_pad[cal_n_hits]/I')
    output_tree.Branch('cal_hit_sector', cal_hit_sector, 'cal_hit_sector[cal_n_hits]/I')
    output_tree.Branch('cal_hit_layer', cal_hit_layer, 'cal_hit_layer[cal_n_hits]/I')
    output_tree.Branch('cal_hit_x', cal_hit_x, 'cal_hit_x[cal_n_hits]/F')
    output_tree.Branch('cal_hit_y', cal_hit_y, 'cal_hit_y[cal_n_hits]/F')
    output_tree.Branch('cal_hit_energy', cal_hit_energy, 'cal_hit_energy[cal_n_hits]/F')

    n_events = input_tree.GetEntries()
    for idx, event in enumerate(input_tree):
        # if idx == 20000:
        #     break

        if idx % (1000) == 0:
            time_min = (time.time() - start_time) // 60
            time_sec = (time.time() - start_time) % 60
            print('Event: {} out of {};'.format(idx, n_events), end=' ')
            print('{} min {} sec'.format(time_min, time_sec))

        hits_tr1, hits_tr2, hits_cal = make_hits_lists(event)

        n_triggers[0] = event.n_triggers
        trigger1[0] = event.trigger1
        trigger2[0] = event.trigger2
        trigger3[0] = event.trigger3

        tr1_n_hits[0] = len(hits_tr1)
        for i, hit in enumerate(hits_tr1):
            tr1_hit_pad[i] = hit.pad
            tr1_hit_sector[i] = hit.sector
            tr1_hit_layer[i] = hit.layer
            tr1_hit_energy[i] = hit.energy
            tr1_hit_x[i] = hit.x
            tr1_hit_y[i] = hit.y

        tr2_n_hits[0] = len(hits_tr2)
        for i, hit in enumerate(hits_tr2):
            tr2_hit_pad[i] = hit.pad
            tr2_hit_sector[i] = hit.sector
            tr2_hit_layer[i] = hit.layer
            tr2_hit_energy[i] = hit.energy
            tr2_hit_x[i] = hit.x
            tr2_hit_y[i] = hit.y

        cal_n_hits[0] = len(hits_cal)
        for i, hit in enumerate(hits_cal):
            cal_hit_pad[i] = hit.pad
            cal_hit_sector[i] = hit.sector
            cal_hit_layer[i] = hit.layer
            cal_hit_x[i] = hit.x
            cal_hit_y[i] = hit.y
            cal_hit_energy[i] = hit.energy

        output_tree.Fill()

    output_tree.Write()
    output_file.Close()

    print("Hooray, extracted tree file is ready, take it :3")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=('Do geometrical selection'))
    parser.add_argument('path_to_file', type=str, help='Provide path to the root file')
    args = parser.parse_args()

    main(args)
