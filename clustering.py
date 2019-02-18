import numpy as np


class Cluster:
    def __init__(self, signals, cluster, weights):
        self.signals = signals
        self.cluster = cluster

        self.energy = self.get_energy()
        self.weights = self.get_weights(weights)

        self.pad = self.get_position('pad')
        self.sector = self.get_position('sector')

        self.rho = self.get_position('rho')

        self.x = self.get_position('x')
        self.y = self.get_position('y')

        self.n_pads = self.get_n_pads()

    def get_energy(self):
        return sum([signal.energy for signal in self.signals if signal.cluster == self.cluster])

    def get_weights(self, weights):
        if weights == 'Energy':
            return [signal.energy for signal in self.signals if signal.cluster == self.cluster]
        elif weights == 'logW':
            w0 = 3.4
            return [max(0, w0 + np.log(signal.energy / self.energy)) for signal in self.signals if signal.cluster == self.cluster]

    def get_position(self, pos_string):
        positions_list = [getattr(signal, pos_string) for signal in self.signals if signal.cluster == self.cluster]
        return sum([positions_list[idx] * weight / sum(self.weights) for idx, weight in enumerate(self.weights)])

    def get_n_pads(self):
        return len([signal for signal in self.signals if signal.cluster == self.cluster])

    def merge(self, cluster2):
            merged_weights = sum(self.weights) + sum(cluster2.weights)
            self.pad = (self.pad * sum(self.weights) + cluster2.pad * sum(cluster2.weights)) / merged_weights
            self.sector = (self.sector * sum(self.weights) + cluster2.sector * sum(cluster2.weights)) / merged_weights
            self.x = (self.x * sum(self.weights) + cluster2.x * sum(cluster2.weights)) / merged_weights
            self.y = (self.y * sum(self.weights) + cluster2.y * sum(cluster2.weights)) / merged_weights
            self.rho = (self.rho * sum(self.weights) + cluster2.rho * sum(cluster2.weights)) / merged_weights
            self.n_pads += cluster2.n_pads
            self.energy += cluster2.energy


def set_neighbors(signals_list):

    for signal in signals_list:
        center_sec, center_pad = signal.sector, signal.pad
        neighbors = []

        for signal_neighbor in signals_list:
            if (signal_neighbor.sector in range(center_sec - 1, center_sec + 2)
               and signal_neighbor.pad in range(center_pad - 1, center_pad + 2)):
                neighbors.append(signal_neighbor)

        neighbors_sorted = sorted(neighbors, key=lambda x: x.energy, reverse=True)
        signal.neighbor = neighbors_sorted[0]


def set_clusters(signals_list):
    cluster_idx = 0
    for signal in signals_list:

        if (signal.neighbor.sector, signal.neighbor.pad) == (signal.sector, signal.pad):
            signal.cluster = cluster_idx
            cluster_idx += 1

    n_non_clusters = -1
    while n_non_clusters != 0:
        n_non_clusters = 0
        for signal in signals_list:
            if signal.cluster == -1:
                n_non_clusters += 1
                if signal.neighbor.cluster != -1:
                    signal.cluster = signal.neighbor.cluster


def merge_clusters(signals_list, cluster_list):

    restart = True
    while restart:
        for cluster1, cluster2 in ((cl1, cl2) for cl1 in cluster_list for cl2 in cluster_list):
            if cluster1 == cluster2:
                continue
            else:
                distance = abs(cluster1.pad - cluster2.pad)
                ratio = cluster2.energy / cluster1.energy

                if distance < 5 or (distance < 20 and ratio < 0.1 - 0.1 / 20 * distance):

                    # This section only to update cluster indices in signals array!
                    # If will be needed to plot 2d map to check clusters
                    cluster1_idx = cluster_list.index(cluster1)
                    cluster2_idx = cluster_list.index(cluster2)
                    for item in signals_list:
                        if item.cluster == cluster2_idx:
                            item.cluster = cluster1_idx
                    for item in signals_list:
                        if item.cluster > cluster2_idx:
                            item.cluster -= 1
                    # End of this section

                    cluster1.merge(cluster2)
                    cluster_list.remove(cluster2)
                    break
        else:
            restart = False
