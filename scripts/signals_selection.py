# File are taken from Sasha Borysov directory:
# /data/alzta/aborysov/tb_2016_data/code/tb16_reco_cd_tr_nn_wfita

#      0           1         2         3
#  _________________________________________
# |   63     |    127   |   191   |   255   |
# |-----------------------------------------|
#  |   62    |    126   |   190   |  254   |
#  |---------------------------------------|
#   |   61    |   125   |   189  |  253   |
#   |-------------------------------------|

#          .........................
#          .........................

#          |  0  |  64 | 128 | 192 |
#          |_____|_____|_____|_____|

# ###CalibFiles###
# Energy calibration files for APVs

# # Important notes which confused me in the past
# # 4 sectors: 0, 1, 2, 3
# # 64 pads: 0, 1, 2, ..., 63
# # 8 layers: 0, 1 - trackers; 2, 3, 4, 5, 6, 7 - calorimeter; 7 - tab (bad)

from ROOT import TFile, TTree, TGraphErrors
import array
import time
import numpy as np
from itertools import islice

import cProfile
import pstats

import argparse

calib_graphs = []


"""Define calibration graphs for APVs"""
# path on alzt.tau.ac.il server = '/data/alzta/aborysov/tb_2016_data/code/lumical_clust/fcalib/'
path = "../apv_calibration/"
for i in range(16):
    calib_file = path + "calibration_apv_{}".format(i) + ".txt"

    # 1st point
    x = [0.]
    y = [0.]
    x_err = [1.e-5]
    y_err = [1.e-5]

    # Calibration x-y data is inverted
    with open(calib_file, 'r') as file:
        for line in islice(file, 1, None):
            x.append(float(line.split('  ')[1]))
            y.append(float(line.split('  ')[0]))
            x_err.append(float(line.split('  ')[3]))
            y_err.append(float(line.split('  ')[2]))

    x = np.array(x)
    # Calibration for tracker APVs are manualy scaled to match MC MPV hit energy
    if i == 0:
        y = np.array(y) * 19.54364863654917
    elif i == 1:
        y = np.array(y) * 18.303542363112417
    elif i == 2:
        y = np.array(y) * 21.093676081159632
    elif i == 3:
        y = np.array(y) * 20.77784418996082
    else:
        y = np.array(y) * 19.206
    x_err = np.array(x_err)
    y_err = np.array(y_err)

    calib_graphs.append(TGraphErrors(len(x), x, y, x_err, y_err))


class ApvMaps:
    '''Define channel to pad conversion'''
    tb15_master = [190 - i if i < 63 else i + 129 for i in range(127)] + [-1]

    tb15_slave = [-1, 62, 63, 60, 61, 58, 59, 56, 57, 54, 55, 52, 53, 50, 51,
                  48, 49, 46, 47, 44, 45, 42, 43, 40, 41, 38, 39, 36, 37, 34,
                  35, 32, 33, 30, 31, 28, 29, 26, 27, 24, 25, 22, 23, 20, 21,
                  18, 19, 16, 17, 14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5,
                  2, 3, 0, 1, 65, 64, 67, 66, 69, 68, 71, 70, 73, 72, 75, 74,
                  77, 76, 79, 78, 81, 80, 83, 82, 85, 84, 87, 86, 89, 88, 91,
                  90, 93, 92, 95, 94, 97, 96, 99, 98, 101, 100, 103, 102, 105,
                  104, 107, 106, 109, 108, 111, 110, 113, 112, 115, 114, 117,
                  116, 119, 118, 121, 120, 123, 122, 125, 124, 127]

    tb15_slave = [-1, 62, 63, 60, 61, 58, 59, 56, 57, 54, 55, 52, 53, 50, 51,
                  48, 49, 46, 47, 44, 45, 42, 43, 40, 41, 38, 39, 36, 37, 34,
                  35, 32, 33, 30, 31, 28, 29, 26, 27, 24, 25, 22, 23, 20, 21,
                  18, 19, 16, 17, 14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2,
                  3, 0, 1, 65, 64, 67, 66, 69, 68, 71, 70, 73, 72, 75, 74, 77,
                  76, 79, 78, 81, 80, 83, 82, 85, 84, 87, 86, 89, 88, 91, 90, 93,
                  92, 95, 94, 97, 96, 99, 98, 101, 100, 103, 102, 105, 104, 107,
                  106, 109, 108, 111, 110, 113, 112, 115, 114, 117, 116, 119, 118,
                  121, 120, 123, 122, 125, 124, 127]

    tb16_master_divider = [-1, 255, 254, 253, 252, 251, 250, 249, 248, 247, 246,
                           245, 244, 243, 242, 241, 240, 239, 238, 237, 236, 235,
                           234, 233, 232, 231, 230, 229, 228, 227, 226, 225, 224,
                           223, 222, 221, 220, 219, 218, 217, 216, 215, 214, 213,
                           212, 211, 210, 209, 208, 207, 206, 205, 204, 203, 202,
                           201, 200, 199, 198, 197, 196, 195, 194, 193, 192, 128,
                           129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139,
                           140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150,
                           151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161,
                           162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172,
                           173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183,
                           184, 185, 186, 187, 188, 189, 191]

    tb16_slave_divider = [126, 124, 125, 122, 123, 120, 121, 118, 119, 116, 117,
                          114, 115, 112, 113, 110, 111, 108, 109, 106, 107, 104,
                          105, 102, 103, 100, 101, 98, 99, 96, 97, 94, 95, 92, 93,
                          90, 91, 88, 89, 86, 87, 84, 85, 82, 83, 80, 81, 78, 79, 76,
                          77, 74, 75, 72, 73, 70, 71, 68, 69, 66, 67, 64, 65, 1, 0, 3,
                          2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14, 17, 16, 19, 18,
                          21, 20, 23, 22, 25, 24, 27, 26, 29, 28, 31, 30, 33, 32, 35,
                          34, 37, 36, 39, 38, 41, 40, 43, 42, 45, 44, 47, 46, 49, 48,
                          51, 50, 53, 52, 55, 54, 57, 56, 59, 58, 61, 60, 63, 62, -1]

    tb16_master_tab_divider = [191, 189, 188, 187, 186, 185, 184, 183, 182, 181, 180,
                               179, 178, 177, 176, 175, 174, 173, 172, 171, 170, 169,
                               168, 167, 166, 165, 164, 163, 162, 161, 160, 159, 158,
                               157, 156, 155, 154, 153, 152, 151, 150, 149, 148, 147,
                               146, 145, 144, 143, 142, 141, 140, 139, 138, 137, 136,
                               135, 134, 133, 132, 131, 130, 129, 128, 192, 193, 194,
                               195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205,
                               206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216,
                               217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227,
                               228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238,
                               239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249,
                               250, 251, 252, 253, 254, 255, -1]

    tb16_slave_tab_divider = [126, 124, 125, 122, 123, 120, 121, 118, 119, 116, 117,
                              114, 115, 112, 113, 110, 111, 108, 109, 106, 107, 104,
                              105, 102, 103, 100, 101, 98, 99, 96, 97, 94, 95, 92, 93,
                              90, 91, 88, 89, 86, 87, 84, 85, 82, 83, 80, 81, 78, 79,
                              76, 77, 74, 75, 72, 73, 70, 71, 68, 69, 66, 67, 64, 65,
                              1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14,
                              17, 16, 19, 18, 21, 20, 23, 22, 25, 24, 27, 26, 29, 28,
                              31, 30, 33, 32, 35, 34, 37, 36, 39, 38, 41, 40, 43, 42,
                              45, 44, 47, 46, 49, 48, 51, 50, 53, 52, 55, 54, 57, 56,
                              59, 58, 61, 60, 63, 62, -1]


def bad_pad(sector, pad, layer):
    """Return true if pad is bad"""
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
            or (layer == 6 and sector == 2 and pad in (34, 42, 54, 57, 59, 60))
            or sector == -1)  # grounded channel


def main(args):
    start_time = time.time()
    # Upload data for analysis
    input_file = TFile.Open(args.path_to_file, "READ")
    input_tree = input_file.apv_reco

    # Create output root file.
    # Create output root file before the tree!!! It prevents memory leakage
    output_file = TFile('./' + args.path_to_file + '_TRANSFORMED.root', "RECREATE")

    # Create output tree
    output_tree = TTree('data', 'Extracted Data')

    n_hits = array.array('i', [0])
    pad = array.array('i', [0] * 64 * 4 * 8)
    sector = array.array('i', [0] * 64 * 4 * 8)
    layer = array.array('i', [0] * 64 * 4 * 8)
    energy = array.array('f', [0.0] * 64 * 4 * 8)

    # Create branches in the output tree for these variables
    output_tree.Branch('n_hits', n_hits, 'n_hits/I')
    output_tree.Branch('pad', pad, 'pad[n_hits]/I')
    output_tree.Branch('sector', sector, 'sector[n_hits]/I')
    output_tree.Branch('layer', layer, 'layer[n_hits]/I')
    output_tree.Branch('energy', energy, 'energy[n_hits]/F')

    n_events = input_tree.GetEntries()
    for idx, event in enumerate(input_tree):
        # if idx == 10000:
        #     break
        if idx % (1000) == 0:
            time_min = (time.time() - start_time) // 60
            time_sec = (time.time() - start_time) % 60
            print('Event: {} out of {};'.format(idx, n_events), end=' ')
            print('{} min {} sec'.format(time_min, time_sec))

        j = 0
        for apv_id, apv_ch, signal, nn, tau, t0, t1 in zip(event.apv_id,
                                                           event.apv_ch,
                                                           event.apv_signal_maxfit,
                                                           event.apv_nn_output,
                                                           event.apv_fit_tau,
                                                           event.apv_fit_t0,
                                                           event.apv_bint1):

            if (tau < 1 or tau > 3 or signal < 0. or signal > 2000. or t0 < (t1 - 2.7) or t0 > (t1 - 0.5) or nn < 0.5):
                continue

            if apv_id < 4:
                apv_map = ApvMaps.tb15_slave if apv_id % 2 == 1 else ApvMaps.tb15_master
            elif 4 <= apv_id < 14:
                apv_map = ApvMaps.tb16_slave_divider if apv_id % 2 == 1 else ApvMaps.tb16_master_divider
            elif apv_id == 14:
                apv_map = ApvMaps.tb16_master_tab_divider
            elif apv_id == 15:
                apv_map = ApvMaps.tb16_slave_tab_divider

            sector[j] = apv_map[apv_ch] // 64
            pad[j] = apv_map[apv_ch] % 64
            layer[j] = apv_id // 2
            if bad_pad(sector[j], pad[j], layer[j]):
                continue
            # Return hit's energy in MIP
            energy[j] = calib_graphs[apv_id].Eval(signal if signal < 1450. else 1450.)
            j += 1

        n_hits[0] = j
        output_tree.Fill()

    output_tree.Write()
    output_file.Close()

    print("Hooray, extracted tree file is ready, take it :3")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=('Signals selection'))
    parser.add_argument('path_to_file', type=str, help='Provide path to the root file')
    args = parser.parse_args()

    pr = cProfile.Profile()
    pr.enable()
    main(args)
    pr.disable()

    ps = pstats.Stats(pr).sort_stats(pstats.SortKey.CUMULATIVE)
    ps.print_stats()
