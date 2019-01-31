
class Cluster:
    '''
    Collected cluster in calorimeter
    Attributes:
    pad, sector - position of the cluster (TOWERS were used for clustering)
    x, y - same, but in Cartesian coorinates
    energy - total energy of all cells in cluster
    n_pads - number of towers which create a cluster

    get_energy() - calculates cluster energy based on signal's data
    get_position() - calculates cluster position based on signal's data
    get_n_pads() - calculates number of pads in cluster based on signal's data
    merge() - update properties of cluster "self" if it was merged with cluster2.
    Sum energy and number of pads. Calculate new weighted average position.
    Do nothing with cluster2. Should be deleted in the code!
    '''
    def __init__(self, data, cluster):
        self.energy = self.get_energy(data, cluster)

        self.pad = self.get_position('pad', data, cluster)
        self.sector = self.get_position('sector', data, cluster)
        self.x = self.get_position('x', data, cluster)
        self.y = self.get_position('y', data, cluster)
        self.n_pads = self.get_n_pads(data, cluster)

    def get_energy(self, data, cluster):
        return sum([signal.energy for signal in data if signal.cluster == cluster])

    def get_position(self, position, data, cluster):
        '''Calculate position as sum with weights(energies) over all points'''
        pos = 0
        pos_energy_list = [(getattr(signal, position), signal.energy) for signal in data if signal.cluster == cluster]
        for pos_energy in pos_energy_list:
            pos += pos_energy[0]*pos_energy[1]/self.energy
        return pos

    def get_n_pads(self, data, cluster):
        return len([signal for signal in data if signal.cluster == cluster])

    def merge(self, cluster2):
            merged_energy = self.energy+cluster2.energy
            self.pad = (self.pad*self.energy+cluster2.pad*cluster2.energy)/merged_energy
            self.sector = (self.sector*self.energy+cluster2.sector*cluster2.energy)/merged_energy
            self.x = (self.x*self.energy+cluster2.x*cluster2.energy)/merged_energy
            self.y = (self.y*self.energy+cluster2.y*cluster2.energy)/merged_energy
            self.n_pads = self.n_pads + cluster2.n_pads
            self.energy = merged_energy


def clustering_in_towers(signals_list,  merge='on'):
    '''
    Group signals into clusters changing their 'cluster' attribute.
    '''
    clusters_list = []

    # Find for every signal the most energetic neighbor
    # And write it as attribute 'neighbor'
    if len(signals_list) == 0:
        return []

    # Sets neighbor for each signal
    set_neighbors(signals_list)

    # Sets cluster index for each signal
    set_clusters(signals_list)

    # Calculate number of primary clusters
    n_clusters = max([signal.cluster for signal in signals_list])+1
    # Add cluster objects to the list.
    for cluster in range(n_clusters):
        clusters_list.append(Cluster(signals_list, cluster))

    # If merge clusters option is on: merge clusters
    if merge == 'on':
        merge_clusters(signals_list, clusters_list)

    # Sort clusters by energy
    clusters_list = sorted(clusters_list, key=lambda x: x.energy, reverse=True)

    # If everything is ok
    return clusters_list


def set_neighbors(signals_list):
    '''
    Finds the most energetic neighbor among neighbors
    '''

    # Loop through all signals for which we are finding the neighbor:
    for signal in signals_list:
        # define position of this signal
        center_sec, center_pad = signal.sector, signal.pad
        # create empty list for it's neighbors
        neighbors = []

        # loop over all neighbor signal-candidate
        for signal_neighbor in signals_list:
            # If signal-candidate is in neighborhood to signal add it to neighbors list.
            # Neighbor can be signal itself! If it is the most energetic among it's neighbors
            if (signal_neighbor.sector in range(center_sec-1, center_sec+2)
               and signal_neighbor.pad in range(center_pad-1, center_pad+2)):
                neighbors.append(signal_neighbor)

        # Sort neighbors by energy
        neighbors_sorted = sorted(neighbors, key=lambda x: x.energy, reverse=True)

        # Pass the highest energetic neighbor to signal attribute "neighbor"
        signal.neighbor = neighbors_sorted[0]


def set_clusters(signals_list):
    '''
    Linking local neighbor algorithm.
    '''

    cluster_idx = 0
    # For signal in data list
    for signal in signals_list:
        # If the neighbor of signal1 is signal1 itself. (This means it is a local maximum)
        # Mark it as a seed (give it cluster_index = 0, 1, 2, ...)
        if (signal.neighbor.sector, signal.neighbor.pad) == (signal.sector, signal.pad):
            signal.cluster = cluster_idx
            cluster_idx += 1

    # Variable to count number of non-cluster sigmals left
    n_non_clusters = -1
    # Stop when all signals got some cluster_index
    while n_non_clusters != 0:
        n_non_clusters = 0
        # for signal in data list
        for signal in signals_list:
            # If signal is not in a cluster: +1 to non_cluster counter
            # If the highest energetic neighbor of signal is in cluster
            # Add this signal to the same cluster
            if signal.cluster == -1:
                n_non_clusters += 1
                if signal.neighbor.cluster != -1:
                    signal.cluster = signal.neighbor.cluster


def merge_clusters(signals_list, cluster_list):
    '''
    Merge clusters if they fulfill 'if' statement condition
    '''

    # restart 'for' loops if clusters merged
    restart = True
    # while clusters merge
    while restart:
        # for cluster1 and cluster2: one line nested double loop
        for cluster1, cluster2 in ((cl1, cl2) for cl1 in cluster_list for cl2 in cluster_list):
            # If it is the same cluster - skip
            if cluster1 == cluster2:
                continue
            else:
                # Calculate distance and energy ratio between clusters
                distance = abs(cluster1.pad - cluster2.pad)
                ratio = cluster2.energy/cluster1.energy

                # Define if statement - when to merge clusters
                if distance < 5 or (distance < 20 and ratio < 0.1-0.1/20*distance):
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

                    # Update cluster1 position, energy, etc.
                    cluster1.merge(cluster2)
                    # Delete 2nd cluster from the list
                    cluster_list.remove(cluster2)
                    # Restart double for loop if clusters merged
                    break
        # If there was no break: no clusters are merged during double for loop.
        # Exit the while loop. All possible clusters already merged
        else:
            restart = False
