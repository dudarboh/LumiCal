
class Detector:
    def __init__(self, signals, clusters):
        self.signals = signals
        self.clusters = clusters

    def n_clusters(self):
        return len(self.clusters)

    def cluster_energy(self, cluster):
        return self.clusters[cluster].energy

    def cluster_pad(self, cluster):
        return self.clusters[cluster].pad

    def cluster_npads(self, cluster):
        return self.clusters[cluster].n_pads

    def cluster_distance(self, cluster1, cluster2):
        return self.clusters[cluster1].pad - self.clusters[cluster2].pad


class Calorimeter(Detector):
    def fill_histos(self, h_dict):
        h_dict['CalNclusters'].Fill(self.n_clusters())
        h_dict['CalClusterEnergy'].Fill(self.cluster_energy(0))
        h_dict['CalClusterPad'].Fill(self.cluster_pad(0))
        h_dict['CalClusterNPads'].Fill(self.cluster_npads(0))


class Tracker(Detector):
    def fill_histos(self, h_dict, tracker_number):
        h_dict[''.join([tracker_number, 'Nclusters'])].Fill(self.n_clusters())
        if self.n_clusters() >= 1:
            h_dict[''.join([tracker_number, 'Cluster1Energy'])].Fill(self.cluster_energy(0))
            h_dict[''.join([tracker_number, 'Cluster1Pad'])].Fill(self.cluster_pad(0))
            h_dict[''.join([tracker_number, 'Cluster1NPads'])].Fill(self.cluster_npads(0))
            if self .n_clusters() >= 2:
                h_dict[''.join([tracker_number, 'Cluster2Energy'])].Fill(self.cluster_energy(1))
                h_dict[''.join([tracker_number, 'Cluster2Pad'])].Fill(self.cluster_pad(1))
                h_dict[''.join([tracker_number, 'Cluster2NPads'])].Fill(self.cluster_npads(1))
                h_dict[''.join([tracker_number, 'Clst1Clst2Distance'])].Fill(self.cluster_distance(1, 0))
        # h_dict[''.join([det_type, 'Cluster1Residuals'])].Fill()
        # h_dict[''.join([det_type, 'Cluster2Residuals'])].Fill()
        # h_dict[''.join([det_type, 'Clst1Clst2ERatio'])].Fill()
        # h_dict[''.join([det_type, 'EResiduals'])].Fill()
        # h_dict[''.join([det_type, 'Clst1Clst2ERatioDistance'])].Fill()
