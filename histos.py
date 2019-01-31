'''
Lists of created histograms in format:
For 1D histograms: (name, title, n_bins, x_min, x_max)
For 2D histograms: (name, title, n_xbins, x_min, x_max, n_ybins, y_min, y_max)
Dont forget to add fill method in the code to fill histograms!
Comment unnecesary!
'''
from ROOT import TH1F, TH2F

histos_1d = [  # N clusters
             ('CalNclusters', 'Calorimeter;N_{clst};N_{ev}', 15, 0, 15),
             ('Tr1Nclusters', 'Tracker1;N_{clst};N_{ev}', 8, 0, 8),
             ('Tr2Nclusters', 'Tracker2;N_{clst};N_{ev}', 8, 0, 8),
             # Calorimeter shower cluster
             ('CalClusterEnergy', 'Shower cluster;E, [MIP];N_{ev}', 2000, 0, 500),
             ('CalClusterPad', 'Shower cluster;Pad;N_{ev}', 64, 0, 64),
             ('CalClusterNPads', 'Shower cluster;N_{pads};N_{ev}', 30, 0, 30),
             # Tracker1 clusters(1,2) properties
             ('Tr1Cluster1Energy', 'Tracker1 cluster1;E, [MIP];N_{ev}', 200, 0, 10),
             ('Tr1Cluster1Pad', 'Tracker1 cluster1;Pad;N_{ev}', 64, 0, 64),
             ('Tr1Cluster1NPads', 'Tracker1 cluster1;N_{pads};N_{ev}', 8, 0, 8),
             ('Tr1Cluster2Energy', 'Tracker1 cluster2;E, [MIP];N_{ev}', 200, 0, 10),
             ('Tr1Cluster2Pad', 'Tracker1 cluster2;Pad;N_{ev}', 64, 0, 64),
             ('Tr1Cluster2NPads', 'Tracker1 cluster2;N_{pads};N_{ev}', 8, 0, 8),
             # Tracker2 clusters(1,2) properties
             ('Tr2Cluster1Energy', 'Tracker2 cluster1;E, [MIP];N_{ev}', 200, 0, 10),
             ('Tr2Cluster1Pad', 'Tracker2 cluster1;Pad;N_{ev}', 64, 0, 64),
             ('Tr2Cluster1NPads', 'Tracker2 cluster1;N_{pads};N_{ev}', 8, 0, 8),
             ('Tr2Cluster2Energy', 'Tracker2 cluster2;E, [MIP];N_{ev}', 200, 0, 10),
             ('Tr2Cluster2Pad', 'Tracker2 cluster2;Pad;N_{ev}', 64, 0, 64),
             ('Tr2Cluster2NPads', 'Tracker2 cluster2;N_{pads};N_{ev}', 8, 0, 8),
             # Residuals of clusters in trackers
             ('Tr1Cluster1Residuals', 'Tracker1 cluster1;d_{clst_fit};N_{ev}', 128, -64, 64),
             ('Tr1Cluster2Residuals', 'Tracker1 cluster2;d_{clst_fit};N_{ev}', 128, -64, 64),
             ('Tr2Cluster1Residuals', 'Tracker2 cluster1;d_{clst_fit};N_{ev}', 128, -64, 64),
             ('Tr2Cluster2Residuals', 'Tracker2 cluster2;d_{clst_fit};N_{ev}', 128, -64, 64),
             # Distance between cluster1 and cluster2 in trackers
             ('Tr1Clst1Clst2Distance', 'Tracker1 clusters:1-2;d_{12};N_{ev}', 128, -64, 64),
             ('Tr2Clst1Clst2Distance', 'Tracker2 clusters:1-2;d_{12};N_{ev}', 128, -64, 64),
             # Energy ratio between cluster1 and cluster2 in trackers
             ('Tr1Clst1Clst2ERatio', 'Tracker1 clusters:1-2;E_{2/1} ratio;N_{ev}', 200, 0, 2),
             ('Tr2Clst1Clst2ERatio', 'Tracker2 clusters:1-2;E_{2/1} ratio;N_{ev}', 200, 0, 2)
            ]

histos_2d = [  # Clusters Energy vs residual in trackers
             ('Tr1EResiduals', 'Tracker1;d_{clst_fit};E, [MIP]', 128, -64, 64, 100, 0, 10),
             ('Tr2EResiduals', 'Tracker2;d_{clst_fit};E, [MIP]', 128, -64, 64, 100, 0, 10),
             # Clusters energy ratio vs distance in trackers
             ('Tr1Clst1Clst2ERatioDistance', 'Tracker1 clusters 1-2;d_{12};E_{2/1} ratio', 128, -64, 64, 200, 0, 2),
             ('Tr2Clst1Clst2ERatioDistance', 'Tracker2 clusters 1-2;d_{12};E_{2/1} ratio', 128, -64, 64, 200, 0, 2)
            ]


def histo_creator(list_1d, list_2d):
    h_dict = {}
    for item in list_1d:
        name, title, n_bins, x_min, x_max = item
        h_dict[name] = TH1F(name, title, n_bins, x_min, x_max)
    for item in list_2d:
        name, title, n_xbins, x_min, x_max, n_ybins, y_min, y_max = item
        h_dict[name] = TH2F(name, title, n_xbins, x_min, x_max, n_ybins, y_min, y_max)
    return h_dict

h_dict = histo_creator(histos_1d, histos_2d)
