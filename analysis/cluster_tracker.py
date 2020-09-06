from ROOT import TFile, TTree
import array
import time
import numpy as np

import argparse


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


def make_clusters_lists(hits_tr1, hits_tr2):
    clusters_tr1 = make_tr_clusters(hits_tr1)
    clusters_tr2 = make_tr_clusters(hits_tr2)
    return clusters_tr1, clusters_tr2


class TrHit:
    def __init__(self, sector, pad, energy):
        self.sector = sector
        self.pad = pad
        self.energy = energy
        self.seed = -1
        rho = 80. + 0.9 + 1.8 * self.pad
        phi = np.pi / 2 + np.pi / 12 - np.pi / 48 - np.pi / 24 * self.sector
        self.x = rho * np.cos(phi)
        self.y = rho * np.sin(phi)


def make_hits_lists(event):
    hits_tracker1 = []
    hits_tracker2 = []

    for (tr1_s, tr1_p, tr1_e, tr2_s, tr2_p, tr2_e) in zip(event.tr1_sector, event.tr1_pad, event.tr1_energy,
                                                        event.tr2_sector, event.tr2_pad, event.tr2_energy):
        hits_tracker1.append(TrHit(tr1_s, tr1_p, tr1_e))
        hits_tracker2.append(TrHit(tr2_s, tr2_p, tr2_e))

        # Selection old
        # if (pad < 20 or sector == 0 or sector == 3
        #     or energy <= 0. or bad_pad(sector, pad, layer)):
        #     continue

    return hits_tracker1, hits_tracker2



def main(args):
    start_time = time.time()

    # Upload data for analysis
    input_file = TFile.Open(args.path_to_file, "READ")
    input_tree = input_file.lumical
    print("Total n events in loaded files: ", input_tree.GetEntries())

    # Create output root file.
    # Create output root file before the tree!!! It prevents memory leakage
    output_file = TFile("./lucas_5gev_e_tr_clusters.root", "RECREATE")
    output_tree = TTree('lumical', 'Clusters in trackers')

    # Create variables associated with the output tree
    tr1_n_clusters = array.array('i', [0])
    tr1_cluster_n_pads = array.array('i', [0] * 128)
    tr1_cluster_pad = array.array('f', [0.0] * 128)
    tr1_cluster_sector = array.array('f', [0.0] * 128)
    tr1_cluster_x = array.array('f', [0.0] * 128)
    tr1_cluster_y = array.array('f', [0.0] * 128)
    tr1_cluster_energy = array.array('f', [0.0] * 128)

    tr2_n_clusters = array.array('i', [0])
    tr2_cluster_n_pads = array.array('i', [0] * 128)
    tr2_cluster_pad = array.array('f', [0.0] * 128)
    tr2_cluster_sector = array.array('f', [0.0] * 128)
    tr2_cluster_x = array.array('f', [0.0] * 128)
    tr2_cluster_y = array.array('f', [0.0] * 128)
    tr2_cluster_energy = array.array('f', [0.0] * 128)

    # Create branches in the output tree for these variables
    output_tree.Branch('tr1_n_clusters', tr1_n_clusters, 'tr1_n_clusters/I')
    output_tree.Branch('tr1_cluster_n_pads', tr1_cluster_n_pads, 'tr1_cluster_n_pads[tr1_n_clusters]/I')
    output_tree.Branch('tr1_cluster_pad', tr1_cluster_pad, 'tr1_cluster_pad[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_sector', tr1_cluster_sector, 'tr1_cluster_sector[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_x', tr1_cluster_x, 'tr1_cluster_x[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_y', tr1_cluster_y, 'tr1_cluster_y[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_energy', tr1_cluster_energy, 'tr1_cluster_energy[tr1_n_clusters]/F')

    output_tree.Branch('tr2_n_clusters', tr2_n_clusters, 'tr2_n_clusters/I')
    output_tree.Branch('tr2_cluster_n_pads', tr2_cluster_n_pads, 'tr2_cluster_n_pads[tr2_n_clusters]/I')
    output_tree.Branch('tr2_cluster_pad', tr2_cluster_pad, 'tr2_cluster_pad[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_sector', tr2_cluster_sector, 'tr2_cluster_sector[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_x', tr2_cluster_x, 'tr2_cluster_x[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_y', tr2_cluster_y, 'tr2_cluster_y[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_energy', tr2_cluster_energy, 'tr2_cluster_energy[tr2_n_clusters]/F')

    n_events = input_tree.GetEntries()
    for idx, event in enumerate(input_tree):
        # if idx == 10000:
        #     break

        if idx % (10000) == 0:
            time_min = (time.time() - start_time) // 60
            time_sec = (time.time() - start_time) % 60
            print('Event: {} out of {};'.format(idx, n_events), end=' ')
            print('{} min {} sec'.format(time_min, time_sec))

        hits_tr1, hits_tr2 = make_hits_lists(event)
        clusters_tr1, clusters_tr2 = make_clusters_lists(hits_tr1, hits_tr2)

        tr1_n_clusters[0] = len(clusters_tr1)
        for i, cluster in enumerate(clusters_tr1):
            tr1_cluster_n_pads[i] = cluster.n_pads
            tr1_cluster_pad[i] = cluster.pad
            tr1_cluster_sector[i] = cluster.sector
            tr1_cluster_x[i] = cluster.x
            tr1_cluster_y[i] = cluster.y
            tr1_cluster_energy[i] = cluster.energy

        tr2_n_clusters[0] = len(clusters_tr2)
        for i, cluster in enumerate(clusters_tr2):
            tr2_cluster_n_pads[i] = cluster.n_pads
            tr2_cluster_pad[i] = cluster.pad
            tr2_cluster_sector[i] = cluster.sector
            tr2_cluster_x[i] = cluster.x
            tr2_cluster_y[i] = cluster.y
            tr2_cluster_energy[i] = cluster.energy


        output_tree.Fill()

    output_tree.Write()
    output_file.Close()

    print("Hooray, extracted tree file is ready, take it :3")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=('Do selection and clustering of data/MC'))
    parser.add_argument('path_to_file', type=str, help='Provide path to the root file')
    args = parser.parse_args()

    main(args)
