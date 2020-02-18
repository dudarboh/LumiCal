import numpy as np
import os.path
from itertools import product
from ROOT import TMath, TCanvas, TH2F, gPad, gStyle, TColor
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



class Tower:
    def __init__(self, tower_hits):
        self.hits = tower_hits
        self.sector = tower_hits[0].sector
        self.pad = tower_hits[0].pad
        self.energy = sum(hit.energy for hit in tower_hits)
        self.n_pads = len(tower_hits)

        self.seed = -1


def make_towers_list(hits_list):
    """Return tower list out of hits list"""
    towers_pos = set((hit.sector, hit.pad) for hit in hits_list)
    towers = []
    for pos in towers_pos:
        tower_hits = [hit for hit in hits_list if (hit.sector, hit.pad) == pos]
        towers.append(Tower(tower_hits))

    towers.sort(key=lambda x: x.energy, reverse=True)

    return towers


def set_tower_seeds(towers_list):
    """Assign cluster indices to the towers. Make towers seeds"""
    seed_idx = 0
    for tower in towers_list:
        if tower.n_pads == 1:
            continue
        sector = tower.sector
        pad = tower.pad

        # Check if it is local maximum
        for tower_neighbor in towers_list:
            if (tower_neighbor.sector in range(sector - 1, sector + 2)
               and tower_neighbor.pad in range(pad - 1, pad + 2)
               and tower_neighbor.energy > tower.energy):
                break
        else:
            # This is local maximum with more than 1 pad.
            tower.seed = seed_idx
            seed_idx += 1


# Assign others towers to seeds
def find_neighbor_assign_cluster(towers_list):
    if not any(tower.seed != -1 for tower in towers_list):
        return
    r = 1
    while any(tower.seed == -1 for tower in towers_list):
        changing = True
        while changing:
            changing = False
            for tower in towers_list:
                # plot_sensor_clusters(towers_list)
                if tower.seed != -1:
                    continue
                sector = tower.sector
                pad = tower.pad

                neighbors = []
                for tower_neighbor in towers_list:
                    if (tower_neighbor.sector in range(sector - r, sector + r + 1)
                       and tower_neighbor.pad in range(pad - r, pad + r + 1)
                       and tower_neighbor.seed != -1):
                        neighbors.append(tower_neighbor)

                # It has neighbors assigned to the cluster
                if len(neighbors) > 0:
                    neighbors.sort(key=lambda x: x.energy, reverse=True)
                    tower.seed = neighbors[0].seed
                    changing = True
        r += 1


class CalCluster:
    def __init__(self, cluster_hits, n_towers):
        self.hits = cluster_hits
        self.energy = sum(hit.energy for hit in cluster_hits)
        self.n_towers = n_towers
        self.n_pads = len(cluster_hits)

        weights = [max(0, 3.4 + np.log(hit.energy / self.energy)) for hit in cluster_hits]
        sum_of_weights = sum(weights)
        if sum_of_weights != 0:
            self.sector = sum(getattr(hit, "sector") * weight for hit, weight in zip(cluster_hits, weights)) / sum_of_weights
            self.pad = sum(getattr(hit, "pad") * weight for hit, weight in zip(cluster_hits, weights)) / sum_of_weights
            self.layer = sum(getattr(hit, "layer") * weight for hit, weight in zip(cluster_hits, weights)) / sum_of_weights
            self.x = sum(getattr(hit, "x") * weight for hit, weight in zip(cluster_hits, weights)) / sum_of_weights
            self.y = sum(getattr(hit, "y") * weight for hit, weight in zip(cluster_hits, weights)) / sum_of_weights
        else:
            self.sector = -999
            self.pad = -999
            self.layer = -999
            self.x = -999
            self.y = -999

    def merge(self, cluster2):
        self.hits.extend(cluster2.hits)
        self.energy = sum(hit.energy for hit in self.hits)
        self.n_towers += cluster2.n_towers
        self.n_pads = len(self.hits)

        weights = [max(0, 3.4 + np.log(hit.energy / self.energy)) for hit in self.hits]
        sum_of_weights = sum(weights)
        if sum_of_weights != 0:
            self.sector = sum(getattr(hit, "sector") * weight for hit, weight in zip(self.hits, weights)) / sum_of_weights
            self.pad = sum(getattr(hit, "pad") * weight for hit, weight in zip(self.hits, weights)) / sum_of_weights
            self.layer = sum(getattr(hit, "layer") * weight for hit, weight in zip(self.hits, weights)) / sum_of_weights
            self.x = sum(getattr(hit, "x") * weight for hit, weight in zip(self.hits, weights)) / sum_of_weights
            self.y = sum(getattr(hit, "y") * weight for hit, weight in zip(self.hits, weights)) / sum_of_weights
        else:
            self.sector = -999
            self.pad = -999
            self.layer = -999
            self.x = -999
            self.y = -999


def merge_clusters(clusters_list):
    """Merge pair of clusters if they meet following condition"""
    while True:
        for clst1, clst2 in product(clusters_list, clusters_list):
            if clst1 == clst2:
                continue
            distance = abs(clst1.y - clst2.y)
            ratio = clst2.energy / clst1.energy
            if distance < 7.5 or (ratio < 0.032 * (20. - distance)):
                clst1.merge(clst2)
                clusters_list.remove(clst2)
                break
        else:
            break


def make_cal_clusters(hits_list):
    """Return cluster list out of towers list"""
    towers_list = make_towers_list(hits_list)

    # plot_sensor_energies(towers_list)

    # This part manages clustering for the calorimeter
    set_tower_seeds(towers_list)

    # plot_sensor_clusters(towers_list)

    find_neighbor_assign_cluster(towers_list)

    clusters = []
    if not towers_list:
        return clusters

    n_clusters = max(tower.seed for tower in towers_list) + 1
    for i in range(n_clusters):
        cluster_hits = []
        n_towers = 0
        for tower in towers_list:
            if tower.seed == i:
                n_towers += 1
                cluster_hits.extend(tower.hits)
        clusters.append(CalCluster(cluster_hits, n_towers))

    clusters.sort(key=lambda x: x.energy, reverse=True)

    merge_clusters(clusters)

    # plot_sensor_clusters(towers_list)
    clusters.sort(key=lambda x: x.energy, reverse=True)

    return clusters


class TrCluster:
    def __init__(self, cluster_hits):
        self.hits = cluster_hits
        self.n_pads = len(cluster_hits)

        weights = [hit.energy for hit in cluster_hits]
        self.energy = sum(weights)
        if self.energy != 0:
            self.sector = sum(getattr(hit, "sector") * weight for hit, weight in zip(cluster_hits, weights)) / self.energy
            self.pad = sum(getattr(hit, "pad") * weight for hit, weight in zip(cluster_hits, weights)) / self.energy
            self.x = sum(getattr(hit, "x") * weight for hit, weight in zip(cluster_hits, weights)) / self.energy
            self.y = sum(getattr(hit, "y") * weight for hit, weight in zip(cluster_hits, weights)) / self.energy
        else:
            self.sector = -999
            self.pad = -999
            self.x = -999
            self.y = -999


def make_tr_clusters(hits_list):
    seed_idx = 0
    for hit in hits_list:
        sector = hit.sector
        pad = hit.pad

        # Check if it is local maximum
        for hit_neighbor in hits_list:
            if (hit_neighbor.sector in range(sector - 1, sector + 2)
               and hit_neighbor.pad in range(pad - 1, pad + 2)
               and hit_neighbor.energy > hit.energy):
                break
        else:
            # This is local maximum.
            hit.seed = seed_idx
            seed_idx += 1

    while any(hit.seed == -1 for hit in hits_list):
        for hit in hits_list:
            if hit.seed != -1:
                continue

            sector = hit.sector
            pad = hit.pad

            neighbors = []
            for hit_neighbor in hits_list:
                if (hit_neighbor.sector in range(sector - 1, sector + 2)
                   and hit_neighbor.pad in range(pad - 1, pad + 2)
                   and hit_neighbor.seed != -1):
                    neighbors.append(hit_neighbor)

            # It has neighbors, check most energetic seed and assign
            if len(neighbors) > 0:
                neighbors.sort(key=lambda x: x.energy, reverse=True)
                hit.seed = neighbors[0].seed

    clusters = []
    if not hits_list:
        return clusters

    n_clusters = max(hit.seed for hit in hits_list) + 1
    for i in range(n_clusters):
        cluster_hits = []
        for hit in hits_list:
            if hit.seed == i:
                cluster_hits.append(hit)
        clusters.append(TrCluster(cluster_hits))

    clusters.sort(key=lambda x: x.energy, reverse=True)

    return clusters


def make_clusters_lists(hits_tr1, hits_tr2, hits_cal):
    clusters_tr1 = make_tr_clusters(hits_tr1)
    clusters_tr2 = make_tr_clusters(hits_tr2)
    clusters_cal = make_cal_clusters(hits_cal)
    return clusters_tr1, clusters_tr2, clusters_cal





















def plot_sensor_energies(towers_list):

    c = TCanvas("c", "c", 800, 1200)
    # gStyle.SetPalette(0)
    c.SetGridx(1)
    c.SetGridy(1)
    gStyle.SetOptStat(0)
    h = TH2F("h", "title", 4, 0., 4., 44, 20., 64.)
    h.GetXaxis().SetTitle("sector number")
    h.GetYaxis().SetTitle("pad number")
    for t in towers_list:
        h.Fill(t.sector, t.pad, t.energy)

    h.SetTitle("")
    h.DrawCopy("colztext")
    # input("wait plot sensor")
    c.Print("./energies.png")
    TColor.InvertPalette()

def plot_sensor_clusters(towers_list):
    pic_number = 0

    c = TCanvas("c", "c", 800, 1200)
    c.SetGridx(1)
    c.SetGridy(1)

    gStyle.SetOptStat(0)
    h = TH2F("h", "title", 4, 0., 4., 44, 20., 64.)
    h.GetXaxis().SetTitle("sector number")
    h.GetYaxis().SetTitle("pad number")
    for t in towers_list:
        if t.seed != -1:
            h.Fill(t.sector, t.pad, t.seed + 1)
        else:
            h.Fill(t.sector, t.pad, t.seed)

    h.SetTitle("")
    h.DrawCopy("colz text")
    # input("wait plot sensor")
    while os.path.exists("./clustering{}.png".format(pic_number)):
        pic_number += 1

    c.Print("./clustering{}.png".format(pic_number))
    pic_number += 1
