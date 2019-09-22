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


from ROOT import TFile, TTree, TChain, TGraphErrors, gROOT
import array
import time
import numpy as np
from itertools import islice

import cProfile
import pstats

# Create output tree
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
tr1_cluster_pad = array.array('f', [0.0] * 128)
tr1_cluster_sector = array.array('f', [0.0] * 128)
tr1_cluster_layer = array.array('f', [0.0] * 128)
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
tr2_cluster_pad = array.array('f', [0.0] * 128)
tr2_cluster_sector = array.array('f', [0.0] * 128)
tr2_cluster_layer = array.array('f', [0.0] * 128)
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
cal_n_towers = array.array('i', [0])
cal_tower_pad = array.array('i', [0] * 128)
cal_tower_sector = array.array('i', [0] * 128)
cal_tower_energy = array.array('f', [0.0] * 128)
cal_tower_cluster = array.array('i', [0] * 128)
cal_n_clusters = array.array('i', [0])
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
output_tree.Branch('tr1_cluster_pad', tr1_cluster_pad, 'tr1_cluster_pad[tr1_n_clusters]/F')
output_tree.Branch('tr1_cluster_sector', tr1_cluster_sector, 'tr1_cluster_sector[tr1_n_clusters]/F')
output_tree.Branch('tr1_cluster_layer', tr1_cluster_layer, 'tr1_cluster_layer[tr1_n_clusters]/F')
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
output_tree.Branch('tr2_cluster_pad', tr2_cluster_pad, 'tr2_cluster_pad[tr2_n_clusters]/F')
output_tree.Branch('tr2_cluster_sector', tr2_cluster_sector, 'tr2_cluster_sector[tr2_n_clusters]/F')
output_tree.Branch('tr2_cluster_layer', tr2_cluster_layer, 'tr2_cluster_layer[tr2_n_clusters]/F')
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
output_tree.Branch('cal_n_towers', cal_n_towers, 'cal_n_towers/I')
output_tree.Branch('cal_tower_pad', cal_tower_pad, 'cal_tower_pad[cal_n_towers]/I')
output_tree.Branch('cal_tower_sector', cal_tower_sector, 'cal_tower_sector[cal_n_towers]/I')
output_tree.Branch('cal_tower_energy', cal_tower_energy, 'cal_tower_energy[cal_n_towers]/F')
output_tree.Branch('cal_tower_cluster', cal_tower_cluster, 'cal_tower_cluster[cal_n_towers]/I')
output_tree.Branch('cal_n_clusters', cal_n_clusters, 'cal_n_clusters/I')
output_tree.Branch('cal_cluster_pad', cal_cluster_pad, 'cal_cluster_pad[cal_n_clusters]/F')
output_tree.Branch('cal_cluster_sector', cal_cluster_sector, 'cal_cluster_sector[cal_n_clusters]/F')
output_tree.Branch('cal_cluster_layer', cal_cluster_layer, 'cal_cluster_layer[cal_n_clusters]/F')
output_tree.Branch('cal_cluster_x', cal_cluster_x, 'cal_cluster_x[cal_n_clusters]/F')
output_tree.Branch('cal_cluster_y', cal_cluster_y, 'cal_cluster_y[cal_n_clusters]/F')
output_tree.Branch('cal_cluster_energy', cal_cluster_energy, 'cal_cluster_energy[cal_n_clusters]/F')


calib_graphs = []


def get_calib_graphs(cls):
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
            or (layer == 6 and sector == 2 and pad in (34, 42, 54, 57, 59, 60)))


def align_detector(hits_tr1, hits_tr2, hits_cal):
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


class Hit:
    def __init__(self, apv_id, apv_channel, apv_signal):
        self.sector, self.pad, self.layer = self.position(apv_id, apv_channel)
        self.energy = self.calib_energy(apv_id, apv_signal)

        rho = 80. + 0.9 + 1.8 * self.pad
        phi = np.pi / 2 + np.pi / 12 - np.pi / 48 - np.pi / 24 * self.sector

        self.x = rho * np.cos(phi)
        self.y = rho * np.sin(phi)

    def position(self, apv_id, apv_channel):
        """Return hit's position: sector, pad, layer"""
        if apv_id < 4:
            apv_map = ApvMaps.tb15_slave if apv_id % 2 == 1 else ApvMaps.tb15_master

        elif 4 <= apv_id < 14:
            apv_map = ApvMaps.tb16_slave_divider if apv_id % 2 == 1 else ApvMaps.tb16_master_divider

        elif apv_id == 14:
            apv_map = ApvMaps.tb16_master_tab_divider

        elif apv_id == 15:
            apv_map = ApvMaps.tb16_slave_tab_divider

        layer = apv_id // 2
        sector = apv_map[apv_channel] // 64
        pad = apv_map[apv_channel] % 64
        return sector, pad, layer

    def calib_energy(self, apv_id, apv_signal):
        """Return hit's energy in MIP"""
        signal = apv_signal if apv_signal < 1450. else 1450.
        return calib_graphs[apv_id].Eval(signal)


class Tower:
    def __init__(self, tower_hits):
        self.hits = tower_hits

        self.sector = self.hits[0].sector
        self.pad = self.hits[0].pad

        self.energy = sum(hit.energy for hit in self.hits)

        self.neighbor = -1
        self.cluster = -1


class Cluster:
    def __init__(self, cluster_hits, det):
        self.hits = cluster_hits
        self.det = det

        self.energy = self.get_energy()

        # Energies in trackers. LogW in Calorimeter
        self.weights = self.get_weights()
        self.sum_of_weights = sum(self.weights)
        if self.sum_of_weights != 0:
            self.sector = self.get_cluster_pos("sector")
            self.pad = self.get_cluster_pos("pad")
            self.layer = self.get_cluster_pos("layer")
            self.x = self.get_cluster_pos("x")
            self.y = self.get_cluster_pos("y")
        else:
            self.sector = -999
            self.pad = -999
            self.layer = -999
            self.x = -999
            self.y = -999

    def get_energy(self):
        """Return total energy of all hits"""
        return sum(hit.energy for hit in self.hits)

    def get_weights(self):
        """Return list of hits' weights"""
        if self.det == 'Cal':
            w0 = 3.4
            return [max(0, w0 + np.log(hit.energy / self.energy)) for hit in self.hits]
        else:
            return [hit.energy for hit in self.hits]

    def get_cluster_pos(self, pos):
        """Return weighted cluster position coordinate: x, y, pad, sector, layer"""
        return sum(getattr(hit, pos) * weight for hit, weight in zip(self.hits, self.weights)) / self.sum_of_weights

    def merge(self, cluster2):
        self.hits.extend(cluster2.hits)
        self.energy = self.get_energy()
        self.weights = self.get_weights()
        self.sum_of_weights = sum(self.weights)
        if self.sum_of_weights != 0:
            self.sector = self.get_cluster_pos("sector")
            self.pad = self.get_cluster_pos("pad")
            self.layer = self.get_cluster_pos("layer")
            self.x = self.get_cluster_pos("x")
            self.y = self.get_cluster_pos("y")
        else:
            self.sector = -999
            self.pad = -999
            self.layer = -999
            self.x = -999
            self.y = -999


def set_most_energetic_neighbors(towers_list):
    """Assign most energetic neighbor to the towers"""
    for tower in towers_list:
        neighbor = tower
        for tower_neighbor in towers_list:
            if (tower_neighbor.sector in range(tower.sector - 1, tower.sector + 2)
               and tower_neighbor.pad in range(tower.pad - 1, tower.pad + 2)
               and tower_neighbor.energy >= neighbor.energy):
                neighbor = tower_neighbor
        tower.neighbor = neighbor


def set_tower_clusters(towers_list):
    """Assign cluster indices to the towers"""
    cluster_idx = 0
    for tower in towers_list:
        if (tower.neighbor.sector, tower.neighbor.pad) == (tower.sector, tower.pad):
            tower.cluster = cluster_idx
            cluster_idx += 1

    while any(tower.cluster == -1 for tower in towers_list):
        for tower in towers_list:
            if tower.cluster == -1 and tower.neighbor.cluster != -1:
                tower.cluster = tower.neighbor.cluster


def merge_clusters(clusters_list):
    """Merge pair of clusters if they meet following condition"""
    restart = True
    while restart:
        for clst1, clst2 in zip(clusters_list, clusters_list):
            if clst1 == clst2:
                continue
            distance = abs(clst1.y - clst2.y)
            ratio = clst2.energy / clst1.energy
            if distance < 9. or (distance < 36. and ratio < 0.1 - 0.1 / 20 * distance):
                clst1.merge(clst2)
                clusters_list.remove(clst2)
                break
        else:
            restart = False


def make_hits_lists(event):
    """Return hits lists from a root file event"""

    hits_calorimeter = []
    hits_tracker1 = []
    hits_tracker2 = []

    for apv_id, apv_ch, signal, nn, tau, t0, t1 in zip(event.apv_id,
                                                       event.apv_ch,
                                                       event.apv_signal_maxfit,
                                                       event.apv_nn_output,
                                                       event.apv_fit_tau,
                                                       event.apv_fit_t0,
                                                       event.apv_bint1):
        if (tau < 1 or tau > 3
           or signal > 2000.
           or t0 < (t1 - 2.7)
           or t0 > (t1 - 0.5)
           or nn < 0.5):
            continue

        hit = Hit(apv_id, apv_ch, signal)
        sector = hit.sector
        pad = hit.pad
        layer = hit.layer

        if (pad < 20
           or (layer > 1 and hit.energy <= 0.)
           or sector == 0 or sector == 3
           or layer == 7
           or (layer <= 1 and signal < 0.)
           or bad_pad(sector, pad, layer)
           or sector < 0):  # THIS IS ESSENTIAL to exclude grounded channel!!!
            continue

        if layer == 0:
            hits_tracker1.append(hit)
        elif layer == 1:
            hits_tracker2.append(hit)
        else:
            hits_calorimeter.append(hit)

    return hits_tracker1, hits_tracker2, hits_calorimeter


def make_towers_list(hits_list):
    """Return tower list out of hits list"""
    towers_pos = set((hit.sector, hit.pad) for hit in hits_list)
    towers = []
    for pos in towers_pos:
        tower_hits = [hit for hit in hits_list if (hit.sector, hit.pad) == pos]
        towers.append(Tower(tower_hits))

    towers.sort(key=lambda x: x.energy, reverse=True)

    return towers


def make_clusters_list(towers_list, det):
    """Return cluster list out of towers list"""
    if not towers_list:
        return []
    clusters = []

    n_clusters = max(tower.cluster for tower in towers_list) + 1
    for i in range(n_clusters):
        cluster_hits = []
        for tower in towers_list:
            if tower.cluster == i:
                cluster_hits.extend(tower.hits)
        clusters.append(Cluster(cluster_hits, det))

    # Sort to start merging the most energetic ones
    clusters.sort(key=lambda x: x.energy, reverse=True)

    merge_clusters(clusters)
    clusters.sort(key=lambda x: x.energy, reverse=True)

    return clusters


def main(root_file):
    start_time = time.time()

    # Upload data for analysis
    input_file = TFile.Open(root_file, "READ")
    input_tree = input_file.apv_reco

    # Create output root file
    output_file = TFile('../extracted_root_files/extracted_data_RENAME.root', "RECREATE")

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
        align_detector(hits_tr1, hits_tr2, hits_cal)

        towers_tr1 = make_towers_list(hits_tr1)
        towers_tr2 = make_towers_list(hits_tr2)
        towers_cal = make_towers_list(hits_cal)

        if towers_tr1:
            set_most_energetic_neighbors(towers_tr1)
            set_tower_clusters(towers_tr1)
        if towers_tr2:
            set_most_energetic_neighbors(towers_tr2)
            set_tower_clusters(towers_tr2)
        if towers_cal:
            set_most_energetic_neighbors(towers_cal)
            set_tower_clusters(towers_cal)

        clusters_tr1 = make_clusters_list(towers_tr1, det='Tr1')
        clusters_tr2 = make_clusters_list(towers_tr2, det='Tr2')
        clusters_cal = make_clusters_list(towers_cal, det='Cal')

        # Resort clusters in trackers by distance to main cluster in calorimeter
        if len(clusters_cal) != 0:
            clusters_tr1.sort(key=lambda x: abs(x.y - clusters_cal[0].y))
            clusters_tr2.sort(key=lambda x: abs(x.y - clusters_cal[0].y))

        output_tree.tr1_n_hits[0] = len(hits_tr1)
        for i, hit in enumerate(hits_tr1):
            output_tree.tr1_hit_pad[i] = hit.pad
            output_tree.tr1_hit_sector[i] = hit.sector
            output_tree.tr1_hit_layer[i] = hit.layer
            output_tree.tr1_hit_x[i] = hit.x
            output_tree.tr1_hit_y[i] = hit.y
            output_tree.tr1_hit_energy[i] = hit.energy
        output_tree.tr1_n_clusters[0] = len(clusters_tr1)
        for i, cluster in enumerate(clusters_tr1):
            output_tree.tr1_cluster_pad[i] = cluster.pad
            output_tree.tr1_cluster_sector[i] = cluster.sector
            output_tree.tr1_cluster_layer[i] = cluster.layer
            output_tree.tr1_cluster_x[i] = cluster.x
            output_tree.tr1_cluster_y[i] = cluster.y
            output_tree.tr1_cluster_energy[i] = cluster.energy
            output_tree.tr1_cluster_n_pads[i] = cluster.n_pads

        output_tree.tr2_n_hits[0] = len(hits_tr2)
        for i, hit in enumerate(hits_tr2):
            output_tree.tr2_hit_pad[i] = hit.pad
            output_tree.tr2_hit_sector[i] = hit.sector
            output_tree.tr2_hit_layer[i] = hit.layer
            output_tree.tr2_hit_x[i] = hit.x
            output_tree.tr2_hit_y[i] = hit.y
            output_tree.tr2_hit_energy[i] = hit.energy
        output_tree.tr2_n_clusters[0] = len(clusters_tr2)
        for i, cluster in enumerate(clusters_tr2):
            output_tree.tr2_cluster_pad[i] = cluster.pad
            output_tree.tr2_cluster_sector[i] = cluster.sector
            output_tree.tr2_cluster_layer[i] = cluster.layer
            output_tree.tr2_cluster_x[i] = cluster.x
            output_tree.tr2_cluster_y[i] = cluster.y
            output_tree.tr2_cluster_energy[i] = cluster.energy
            output_tree.tr2_cluster_n_pads[i] = cluster.n_pads

        output_tree.cal_n_hits[0] = len(hits_cal)
        for i, hit in enumerate(hits_cal):
            output_tree.cal_hit_pad[i] = hit.pad
            output_tree.cal_hit_sector[i] = hit.sector
            output_tree.cal_hit_layer[i] = hit.layer
            output_tree.cal_hit_x[i] = hit.x
            output_tree.cal_hit_y[i] = hit.y
            output_tree.cal_hit_energy[i] = hit.energy
        output_tree.cal_n_towers[0] = len(towers_cal)
        for i, tower in enumerate(towers_cal):
            output_tree.cal_tower_pad[i] = tower.pad
            output_tree.cal_tower_sector[i] = tower.sector
            output_tree.cal_tower_energy[i] = tower.energy
            output_tree.cal_tower_cluster[i] = tower.cluster
        output_tree.cal_n_clusters[0] = len(clusters_cal)
        for i, cluster in enumerate(clusters_cal):
            output_tree.cal_cluster_pad[i] = cluster.pad
            output_tree.cal_cluster_sector[i] = cluster.sector
            output_tree.cal_cluster_layer[i] = cluster.layer
            output_tree.cal_cluster_x[i] = cluster.x
            output_tree.cal_cluster_y[i] = cluster.y
            output_tree.cal_cluster_energy[i] = cluster.energy
            output_tree.cal_cluster_n_pads[i] = cluster.n_pads

        output_tree.Fill()

    output_tree.Write()
    output_file.Close()

    print("Hooray, extracted tree file is ready, take it :3")


pr = cProfile.Profile()
pr.enable()
# main("../data/run737_741_5gev.root")
# main("../data/run742_744_4gev.root")
# main("../data/run745_746_3gev.root")
# main("../data/run747_748_2gev.root")
main("../data/run749_750_1gev.root")
pr.disable()

ps = pstats.Stats(pr).sort_stats(pstats.SortKey.CUMULATIVE)
ps.print_stats()
