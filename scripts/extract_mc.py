from ROOT import TFile, TTree
import array

import math
import random
import numpy as np

import time
import cProfile
import pstats

# Important notes which confused me in the past
# 4 sectors: 0, 1, 2, 3
# 64 pads: 0, 1, 2, ..., 63
# 8 layers: 0, 1 - trackers; 2, 3, 4, 5, 6, 7 - calorimeter; 7 - tab (bad)


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


def align_detector(hits_tr1, hits_tr2, hits_cal):
    # Was checked for Itamars misalignment.
    tr1_shift = -0.2480096316087952
    tr2_shift = 0.9012092661162399
    cal_shift = -0.6531996345075299

    for hit in hits_tr1:
        hit.y -= tr1_shift
    for hit in hits_tr2:
        hit.y -= tr2_shift
    for hit in hits_cal:
        hit.y -= cal_shift


class Hit:
    def __init__(self, sector, pad, layer, bs, energy_in_mip):
        self.sector = sector
        self.pad = pad
        self.layer = layer
        self.energy = energy_in_mip
        self.bs = bs
        rho = 80. + 0.9 + 1.8 * self.pad
        phi = np.pi / 2 + np.pi / 12 - np.pi / 48 - np.pi / 24 * self.sector
        self.x = rho * np.cos(phi)
        self.y = rho * np.sin(phi)


class Tower:
    def __init__(self, tower_hits):
        self.hits = tower_hits

        self.sector = self.hits[0].sector
        self.pad = self.hits[0].pad

        self.energy = sum([hit.energy for hit in self.hits])

        self.neighbor = -1
        self.cluster = -1


class Cluster:
    def __init__(self, cluster_hits, cluster_towers, det):
        self.hits = cluster_hits
        self.towers = cluster_towers
        self.det = det

        self.energy = self.get_energy()

        # Energies in trackers. LogW in Calorimeter
        self.weights = self.get_weights()

        self.sector = self.get_cluster_sector()
        self.pad = self.get_cluster_pad()
        self.layer = self.get_cluster_layer()
        self.x = self.get_cluster_x()
        self.y = self.get_cluster_y()

        self.n_pads = self.get_n_pads()

    def get_energy(self):
        return sum([hit.energy for hit in self.hits])

    def get_weights(self):
        if self.det == 'Cal':
            w0 = 3.4
            return [max(0, w0 + np.log(hit.energy / self.energy)) for hit in self.hits]
        else:
            return [hit.energy for hit in self.hits]

    def get_cluster_sector(self):
        return sum([self.hits[i].sector * self.weights[i] for i in range(len(self.hits))]) / sum(self.weights)

    def get_cluster_pad(self):
        return sum([self.hits[i].pad * self.weights[i] for i in range(len(self.hits))]) / sum(self.weights)

    def get_cluster_layer(self):
        return sum([self.hits[i].layer * self.weights[i] for i in range(len(self.hits))]) / sum(self.weights)

    def get_cluster_x(self):
        return sum([self.hits[i].x * self.weights[i] for i in range(len(self.hits))]) / sum(self.weights)

    def get_cluster_y(self):
        return sum([self.hits[i].y * self.weights[i] for i in range(len(self.hits))]) / sum(self.weights)

    def get_n_pads(self):
        return len(self.hits)

    def merge(self, cluster2):
        self.hits.extend(cluster2.hits)
        self.towers.extend(cluster2.towers)
        self.energy = self.get_energy()
        self.weights = self.get_weights()

        self.pad = self.get_cluster_pad()
        self.sector = self.get_cluster_sector()
        self.layer = self.get_cluster_layer()
        self.x = self.get_cluster_x()
        self.y = self.get_cluster_y()

        self.n_pads = self.get_n_pads()


class OutputTree:
    def __init__(self):
        self.tree = TTree('mc', 'Extracted MC')

    def define_arrays(self):
        self.tr1_n_hits = array.array('i', [0])
        self.tr1_hit_pad = array.array('i', [0] * 128)
        self.tr1_hit_sector = array.array('i', [0] * 128)
        self.tr1_hit_layer = array.array('i', [0] * 128)
        self.tr1_hit_x = array.array('f', [0.0] * 128)
        self.tr1_hit_y = array.array('f', [0.0] * 128)
        self.tr1_hit_bs = array.array('i', [0] * 128)
        self.tr1_hit_energy = array.array('f', [0.0] * 128)
        self.tr1_n_clusters = array.array('i', [0])
        self.tr1_cluster_pad = array.array('f', [0.0] * 128)
        self.tr1_cluster_sector = array.array('f', [0.0] * 128)
        self.tr1_cluster_layer = array.array('f', [0.0] * 128)
        self.tr1_cluster_x = array.array('f', [0.0] * 128)
        self.tr1_cluster_y = array.array('f', [0.0] * 128)
        self.tr1_cluster_energy = array.array('f', [0.0] * 128)
        self.tr1_cluster_n_pads = array.array('i', [0] * 128)

        self.tr2_n_hits = array.array('i', [0])
        self.tr2_hit_pad = array.array('i', [0] * 128)
        self.tr2_hit_sector = array.array('i', [0] * 128)
        self.tr2_hit_layer = array.array('i', [0] * 128)
        self.tr2_hit_x = array.array('f', [0.0] * 128)
        self.tr2_hit_y = array.array('f', [0.0] * 128)
        self.tr2_hit_bs = array.array('i', [0] * 128)
        self.tr2_hit_energy = array.array('f', [0.0] * 128)
        self.tr2_n_clusters = array.array('i', [0])
        self.tr2_cluster_pad = array.array('f', [0.0] * 128)
        self.tr2_cluster_sector = array.array('f', [0.0] * 128)
        self.tr2_cluster_layer = array.array('f', [0.0] * 128)
        self.tr2_cluster_x = array.array('f', [0.0] * 128)
        self.tr2_cluster_y = array.array('f', [0.0] * 128)
        self.tr2_cluster_energy = array.array('f', [0.0] * 128)
        self.tr2_cluster_n_pads = array.array('i', [0] * 128)

        self.cal_n_hits = array.array('i', [0])
        self.cal_hit_pad = array.array('i', [0] * 128 * 5)
        self.cal_hit_sector = array.array('i', [0] * 128 * 5)
        self.cal_hit_layer = array.array('i', [0] * 128 * 5)
        self.cal_hit_x = array.array('f', [0.0] * 128 * 5)
        self.cal_hit_y = array.array('f', [0.0] * 128 * 5)
        self.cal_hit_bs = array.array('i', [0] * 128)
        self.cal_hit_energy = array.array('f', [0.0] * 128 * 5)
        self.cal_n_towers = array.array('i', [0])
        self.cal_tower_pad = array.array('i', [0] * 128)
        self.cal_tower_sector = array.array('i', [0] * 128)
        self.cal_tower_energy = array.array('f', [0.0] * 128)
        self.cal_tower_cluster = array.array('i', [0] * 128)
        self.cal_n_clusters = array.array('i', [0])
        self.cal_cluster_pad = array.array('f', [0.0] * 128 * 5)
        self.cal_cluster_sector = array.array('f', [0.0] * 128 * 5)
        self.cal_cluster_layer = array.array('f', [0.0] * 128 * 5)
        self.cal_cluster_x = array.array('f', [0.0] * 128 * 5)
        self.cal_cluster_y = array.array('f', [0.0] * 128 * 5)
        self.cal_cluster_energy = array.array('f', [0.0] * 128 * 5)
        self.cal_cluster_n_pads = array.array('i', [0] * 128 * 5)

    def define_branches(self):
        self.tree.Branch('tr1_n_hits', self.tr1_n_hits, 'tr1_n_hits/I')
        self.tree.Branch('tr1_hit_pad', self.tr1_hit_pad, 'tr1_hit_pad[tr1_n_hits]/I')
        self.tree.Branch('tr1_hit_sector', self.tr1_hit_sector, 'tr1_hit_sector[tr1_n_hits]/I')
        self.tree.Branch('tr1_hit_layer', self.tr1_hit_layer, 'tr1_hit_layer[tr1_n_hits]/I')
        self.tree.Branch('tr1_hit_x', self.tr1_hit_x, 'tr1_hit_x[tr1_n_hits]/F')
        self.tree.Branch('tr1_hit_y', self.tr1_hit_y, 'tr1_hit_y[tr1_n_hits]/F')
        self.tree.Branch('tr1_hit_bs', self.tr1_hit_bs, 'tr1_hit_bs[tr1_n_hits]/I')
        self.tree.Branch('tr1_hit_energy', self.tr1_hit_energy, 'tr1_hit_energy[tr1_n_hits]/F')
        self.tree.Branch('tr1_n_clusters', self.tr1_n_clusters, 'tr1_n_clusters/I')
        self.tree.Branch('tr1_cluster_pad', self.tr1_cluster_pad, 'tr1_cluster_pad[tr1_n_clusters]/F')
        self.tree.Branch('tr1_cluster_sector', self.tr1_cluster_sector, 'tr1_cluster_sector[tr1_n_clusters]/F')
        self.tree.Branch('tr1_cluster_layer', self.tr1_cluster_layer, 'tr1_cluster_layer[tr1_n_clusters]/F')
        self.tree.Branch('tr1_cluster_x', self.tr1_cluster_x, 'tr1_cluster_x[tr1_n_clusters]/F')
        self.tree.Branch('tr1_cluster_y', self.tr1_cluster_y, 'tr1_cluster_y[tr1_n_clusters]/F')
        self.tree.Branch('tr1_cluster_energy', self.tr1_cluster_energy, 'tr1_cluster_energy[tr1_n_clusters]/F')
        self.tree.Branch('tr1_cluster_n_pads', self.tr1_cluster_n_pads, 'tr1_cluster_n_pads[tr1_n_clusters]/I')

        self.tree.Branch('tr2_n_hits', self.tr2_n_hits, 'tr2_n_hits/I')
        self.tree.Branch('tr2_hit_pad', self.tr2_hit_pad, 'tr2_hit_pad[tr2_n_hits]/I')
        self.tree.Branch('tr2_hit_sector', self.tr2_hit_sector, 'tr2_hit_sector[tr2_n_hits]/I')
        self.tree.Branch('tr2_hit_layer', self.tr2_hit_layer, 'tr2_hit_layer[tr2_n_hits]/I')
        self.tree.Branch('tr2_hit_x', self.tr2_hit_x, 'tr2_hit_x[tr2_n_hits]/F')
        self.tree.Branch('tr2_hit_y', self.tr2_hit_y, 'tr2_hit_y[tr2_n_hits]/F')
        self.tree.Branch('tr2_hit_bs', self.tr2_hit_bs, 'tr2_hit_bs[tr2_n_hits]/I')
        self.tree.Branch('tr2_hit_energy', self.tr2_hit_energy, 'tr2_hit_energy[tr2_n_hits]/F')
        self.tree.Branch('tr2_n_clusters', self.tr2_n_clusters, 'tr2_n_clusters/I')
        self.tree.Branch('tr2_cluster_pad', self.tr2_cluster_pad, 'tr2_cluster_pad[tr2_n_clusters]/F')
        self.tree.Branch('tr2_cluster_sector', self.tr2_cluster_sector, 'tr2_cluster_sector[tr2_n_clusters]/F')
        self.tree.Branch('tr2_cluster_layer', self.tr2_cluster_layer, 'tr2_cluster_layer[tr2_n_clusters]/F')
        self.tree.Branch('tr2_cluster_x', self.tr2_cluster_x, 'tr2_cluster_x[tr2_n_clusters]/F')
        self.tree.Branch('tr2_cluster_y', self.tr2_cluster_y, 'tr2_cluster_y[tr2_n_clusters]/F')
        self.tree.Branch('tr2_cluster_energy', self.tr2_cluster_energy, 'tr2_cluster_energy[tr2_n_clusters]/F')
        self.tree.Branch('tr2_cluster_n_pads', self.tr2_cluster_n_pads, 'tr2_cluster_n_pads[tr2_n_clusters]/I')

        self.tree.Branch('cal_n_hits', self.cal_n_hits, 'cal_n_hits/I')
        self.tree.Branch('cal_hit_pad', self.cal_hit_pad, 'cal_hit_pad[cal_n_hits]/I')
        self.tree.Branch('cal_hit_sector', self.cal_hit_sector, 'cal_hit_sector[cal_n_hits]/I')
        self.tree.Branch('cal_hit_layer', self.cal_hit_layer, 'cal_hit_layer[cal_n_hits]/I')
        self.tree.Branch('cal_hit_x', self.cal_hit_x, 'cal_hit_x[cal_n_hits]/F')
        self.tree.Branch('cal_hit_y', self.cal_hit_y, 'cal_hit_y[cal_n_hits]/F')
        self.tree.Branch('cal_hit_bs', self.cal_hit_bs, 'cal_hit_bs[cal_n_hits]/I')
        self.tree.Branch('cal_hit_energy', self.cal_hit_energy, 'cal_hit_energy[cal_n_hits]/F')
        self.tree.Branch('cal_n_towers', self.cal_n_towers, 'cal_n_towers/I')
        self.tree.Branch('cal_tower_pad', self.cal_tower_pad, 'cal_tower_pad[cal_n_towers]/I')
        self.tree.Branch('cal_tower_sector', self.cal_tower_sector, 'cal_tower_sector[cal_n_towers]/I')
        self.tree.Branch('cal_tower_energy', self.cal_tower_energy, 'cal_tower_energy[cal_n_towers]/F')
        self.tree.Branch('cal_tower_cluster', self.cal_tower_cluster, 'cal_tower_cluster[cal_n_towers]/I')
        self.tree.Branch('cal_n_clusters', self.cal_n_clusters, 'cal_n_clusters/I')
        self.tree.Branch('cal_cluster_pad', self.cal_cluster_pad, 'cal_cluster_pad[cal_n_clusters]/F')
        self.tree.Branch('cal_cluster_sector', self.cal_cluster_sector, 'cal_cluster_sector[cal_n_clusters]/F')
        self.tree.Branch('cal_cluster_layer', self.cal_cluster_layer, 'cal_cluster_layer[cal_n_clusters]/F')
        self.tree.Branch('cal_cluster_x', self.cal_cluster_x, 'cal_cluster_x[cal_n_clusters]/F')
        self.tree.Branch('cal_cluster_y', self.cal_cluster_y, 'cal_cluster_y[cal_n_clusters]/F')
        self.tree.Branch('cal_cluster_energy', self.cal_cluster_energy, 'cal_cluster_energy[cal_n_clusters]/F')
        self.tree.Branch('cal_cluster_n_pads', self.cal_cluster_n_pads, 'cal_cluster_n_pads[cal_n_clusters]/I')


def set_most_energetic_neighbors(towers_list):
    for tower in towers_list:
        neighbor = tower
        for tower_neighbor in towers_list:
            if (tower_neighbor.sector in range(tower.sector - 1, tower.sector + 2)
               and tower_neighbor.pad in range(tower.pad - 1, tower.pad + 2)):
                neighbor = tower_neighbor if tower_neighbor.energy >= neighbor.energy else neighbor
        tower.neighbor = neighbor


def set_tower_clusters(towers_list):
    cluster_idx = 0
    for tower in towers_list:
        if (tower.neighbor.sector, tower.neighbor.pad) == (tower.sector, tower.pad):
            tower.cluster = cluster_idx
            cluster_idx += 1

    n_non_clusters = -1
    while n_non_clusters != 0:
        n_non_clusters = 0
        for tower in towers_list:
            if tower.cluster == -1:
                n_non_clusters += 1
                if tower.neighbor.cluster != -1:
                    tower.cluster = tower.neighbor.cluster


def merge_clusters(clusters_list):
    restart = True
    while restart:
        for idx1, cluster1 in enumerate(clusters_list):
            for idx2, cluster2 in enumerate(clusters_list):
                if cluster1 != cluster2:
                    distance = abs(cluster1.y - cluster2.y)
                    ratio = cluster2.energy / cluster1.energy
                    if distance < 9. or (distance < 36. and ratio < 0.1 - 0.1 / 20 * distance):
                        cluster1.merge(cluster2)
                        clusters_list.remove(cluster2)
                        break
            else:
                continue
            break
        else:
            restart = False


def make_hits_lists(event):
    n_hits = event.hit_n
    hit_sector = event.hit_sector
    hit_pad = event.hit_pad
    hit_layer = event.hit_layer

    hit_energy_mev = event.hit_energy
    hit_bs = event.hit_bs

    S0 = 0.819
    p1 = 2.166
    p0 = 0.999 / 2.

    hits_tracker1 = []
    hits_tracker2 = []
    hits_calorimeter = []

    for i in range(n_hits):
        hit_energy = hit_energy_mev[i] / 0.0885

        # Implement noise in progress. Need to put not random number!
        # hit_energy = random.gauss(hit_energy, 0.52353509)

        # Selection as in data
        if (hit_pad[i] < 20
            or (hit_layer[i] > 1 and hit_energy < 1.4)
            or hit_layer[i] == 7
            or hit_sector[i] == 0 or hit_sector[i] == 3
            or bad_pad(hit_sector[i], hit_pad[i], hit_layer[i])
            or (hit_layer[i] <= 1 and hit_energy < 0.)):
            continue

        # Calorimeter efficiency simulation
        if hit_layer[i] > 1:
            if random.random() > (1. + math.erf((hit_energy - S0) / p1)) * p0:
                continue

        # If passed the selection create hit and add to corresponding list
        hit = Hit(hit_sector[i], hit_pad[i], hit_layer[i], hit_bs[i], hit_energy)

        if hit.layer == 0:
            hits_tracker1.append(hit)
        elif hit.layer == 1:
            hits_tracker2.append(hit)
        else:
            hits_calorimeter.append(hit)

    return hits_tracker1, hits_tracker2, hits_calorimeter


def make_towers_list(hits_list):
    towers_pos = set([(hit.sector, hit.pad) for hit in hits_list])
    towers = []
    for pos in towers_pos:
        tower_hits = [hit for hit in hits_list if (hit.sector, hit.pad) == pos]
        towers.append(Tower(tower_hits))

    towers.sort(key=lambda x: x.energy, reverse=True)

    return towers


def make_clusters_list(towers_list, det):
    if not towers_list:
        return []

    clusters = []

    n_clusters = max([tower.cluster for tower in towers_list]) + 1
    for clst_idx in range(n_clusters):
        cluster_hits = []
        cluster_towers = []
        for tower in towers_list:
            if tower.cluster == clst_idx:
                cluster_hits.extend(tower.hits)
                cluster_towers.append(tower)
        clusters.append(Cluster(cluster_hits, cluster_towers, det))

    # Sort to start merging the most energetic ones
    clusters.sort(key=lambda x: x.energy, reverse=True)

    merge_clusters(clusters)
    clusters.sort(key=lambda x: x.energy, reverse=True)

    # Change tower clusters indices after resorting by the energy
    for i, cluster in enumerate(clusters):
        for tower in cluster.towers:
            tower.cluster = i

    return clusters


def main():
    start_time = time.time()

    file = TFile.Open('../mc_root_files/lucas_tb16_5gev.root')
    tree = file.LumiCal
    print("Total n events in loaded files: ", tree.GetEntries())

    output_file = TFile("../extracted_root_files/extracted_mc_5gev.root", "RECREATE")
    output_file.cd()

    output_tree = OutputTree()
    output_tree.define_arrays()
    output_tree.define_branches()

    for idx, event in enumerate(tree):
        # if idx == 200000:
        #     break

        if idx % (10000) == 0:
            time_min = (time.time() - start_time) // 60
            time_sec = (time.time() - start_time) % 60
            print('Event: {} out of {};'.format(idx, tree.GetEntries()), end=' ')
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
            output_tree.tr1_hit_bs[i] = hit.bs
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
            output_tree.tr2_hit_bs[i] = hit.bs
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
            output_tree.cal_hit_bs[i] = hit.bs
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

        output_tree.tree.Fill()

    output_tree.tree.Write()
    output_file.Close()

    print("Hooray, extracted tree file is ready, take it :3")


pr = cProfile.Profile()
pr.enable()
main()
pr.disable()

ps = pstats.Stats(pr).sort_stats(pstats.SortKey.CUMULATIVE)
ps.print_stats()
