'''
This script was written to analise data of the LumiCal test beam runs of 2016.
It's input - root files with signals and apvs' information provided by
Sasha's Borisov readout system.
Script has been tested with: root v6.14/06 and python 3.7.2.
Time needed is ~20 minutes for 55000 events. (AT MY OLD LAPTOP!)
Script can do the following:
1)Show 2d map of signals distribution in the event
2)Do clustering of signals in calorimeter using "linking neighbors" algorithm.
3)Calculate and plot various distributions: Number of clusters/Energy/Position/
Distance between clusters/etc. distributions
4) Show 2d map of signals in the both trackers plane.
5) Plot energy distribution in the trackets and fit it with Landau-Gaus convolution function. (BROKEN IN THE MOMENT)
'''

from ROOT import TH1F, TH2F, TF1, TGraphErrors, TFile, gROOT, TCanvas, TColor, gStyle
import time
import numpy as np

from variables import n_sectors, n_pads, apv_maps, calib_path, \
    calib_file_names, langaufun

# To measure time of execution
start_time = time.time()

# Dont draw pictures in the end
gROOT.SetBatch(1)


class Signal():
    '''
    Signal in the detector
    Attributes:
    sector, pad, layer - position of the signal
    energy - energy of the signal
    neighbor - the most energetic neighbor of the signal or itself.
    x, y - coordinates of signal in x, y.
    cluster - cluster index to which the signal was asigned

    position() - method to calculate position of the signal
    based on apv's id and channel where signal was recorded
    calib_energy - method to calculate calibrated energy of the
    signal.It uses: calibration files for this apv, recorded signal
    value of ADC (Voltage) in apv.
    '''

    def __init__(self, sector, pad, energy):
        # Calculate all parameters of the signal when object created
        self.sector = sector
        self.pad = pad
        self.energy = energy

        self.neighbor = -1
        self.cluster = -1
        # Transform pad/sector numbers into x,y coordinates
        rho = 80+0.9+1.8*self.pad
        phi = np.pi/2+np.pi/12-np.pi/48-np.pi/24*self.sector
        self.x = rho*np.cos(phi)
        self.y = rho*np.sin(phi)


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


class EventData:
    '''
    Main analysis class. It reads signals data from input root file.
    Collect clusters in calorimeter from signals using linking-neighbor
    algorithm. Merge created clusters if they satisfy arbitraty chosen condition.
    Fill all created histograms with corresponded data.

    Attributes:
    create_histos() - create histograms to output
    desired distributions
    fill_histos() - fill those histos with parameters of event.
    extract_data() - reads signals from the root file
    clustering_in_towers() - finds the most energetic neighbors for all signals
    (set_neighbors()). Collect signals into clusters (set_clusters()). Merge
    them if(something) (merge_clusters()).
    PlotCheck - obsolete. To check 2d map of the event
    PlotCluterCheck - more obsolete
    '''
    def __init__(self, filename):
        # Create root file and open "filename" file
        self.file = TFile.Open(filename)

        # Read a TTree from this file
        self.tree = self.file.apv_reco

        # Create all needed histograms
        self.create_histos()

    def create_histos(self):
        '''
        Creates histograms needed for the analysis and writes them in dictionary
        with corresponded keys. Feel free to comment unnecessary histograms to
        not create mess in the output file root.
        '''
        self.h_dict = {}
        # Number of clusters
        h_key = 'h_n_clusters'
        h_title = 'N clusters;N Clusters;N Events'
        self.h_dict[h_key] = TH1F(h_key, h_title, 15, 0, 15)

        # Energy of all 1pad clusters in calorimeter
        h_key = 'h_energy_1pad_clusters'
        h_title = 'Energy: 1 pad clusters;Energy [MIP];N events'
        #self.h_dict[h_key] = TH1F(h_key, h_title, 2000, 0, 20)

        # Energy of all hits in tracker1
        h_key = 'h_tracker1_energy'
        h_title = 'Tracker1 Energy;Energy [MIP];N events'
        self.h_dict[h_key] = TH1F(h_key, h_title, 200, 0, 10)

        # Energy of all hits in tracker2
        h_key = 'h_tracker2_energy'
        h_title = 'Tracker2 Energy;Energy [MIP];N events'
        self.h_dict[h_key] = TH1F(h_key, h_title, 200, 0, 10)

        # 2D:Hits weighted with energy in tracker1
        h_key = '2d_map_tracker1_energy'
        h_title = 'energy weighed cells in tracker1;sector;pad'
        self.h_dict[h_key] = TH2F(h_key, h_title, n_sectors+2, 0, n_sectors+2, 64, 0, 64)

        # 2D:Hits weighted with energy in tracker2
        h_key = '2d_map_tracker2_energy'
        h_title = 'energy weighted cells in tracker2;sector;pad'
        self.h_dict[h_key] = TH2F(h_key, h_title, n_sectors+2, 0, n_sectors+2, 64, 0, 64)

        # 2D:Hits in tracker1
        h_key = '2d_map_tracker1_hits'
        h_title = 'hits in tracker1;sector;pad'
        self.h_dict[h_key] = TH2F(h_key, h_title, n_sectors+2, 0, n_sectors+2, 64, 0, 64)

        # 2D:Hits in tracker2
        h_key = '2d_map_tracker2_hits'
        h_title = 'hits in tracker2;sector;pad'
        self.h_dict[h_key] = TH2F(h_key, h_title, n_sectors+2, 0, n_sectors+2, 64, 0, 64)

        # Distance between cluster and hit in tracker1 weighted with hit energy
        h_key = 'h_cluster_tracker1_distance'
        h_title = 'Tracker1 hits vs cluster;Distance cluster-hit [pad];Summed hit energies'
        self.h_dict[h_key] = TH1F(h_key, h_title, 300, 0, 60)

        # Distance between cluster and hit in tracker2 weighted with hit energy
        h_key = 'h_cluster_tracker2_distance'
        h_title = 'Tracker2 hits vs cluster;Distance cluster-hit [pad];Summed hit energies'
        self.h_dict[h_key] = TH1F(h_key, h_title, 300, 0, 60)

        # Hit position for tracker1
        h_key = 'h_cluster_tracker1_position'
        h_title = 'Tracker1 hits position;position [pad];N hits'
        self.h_dict[h_key] = TH1F(h_key, h_title, 64, 0, 64)

        # Hit position for tracker2
        h_key = 'h_cluster_tracker2_position'
        h_title = 'Tracker2 hits position;position [pad];N hits'
        self.h_dict[h_key] = TH1F(h_key, h_title, 64, 0, 64)

        for cluster in range(1):
            # Cluster energy
            h_key = 'h_cluster_energy_{}'.format(cluster+1)
            h_title = 'Energy: {} clust;Energy [MIP];N events'.format(cluster+1)
            self.h_dict[h_key] = TH1F(h_key, h_title, 2000, 0, 500)

            # Cluster pad position
            h_key = 'h_cluster_pad_pos_{}'.format(cluster+1)
            h_title = 'Position: {} cluster;pos [pad];N events'.format(cluster+1)
            self.h_dict[h_key] = TH1F(h_key, h_title, 200, 0, n_pads)

            # Cluster number of pads
            h_key = 'h_cluster_npads_{}'.format(cluster+1)
            h_title = 'N pads: {} cluster;N pads;N events'.format(cluster+1)
            self.h_dict[h_key] = TH1F(h_key, h_title, 30, 0, 30)

        # for cluster1 in range(1):
        #     for cluster2 in range(1, 3):
        #         # Distance in pads between cluster1 and cluster2
        #         h_key = 'h_cluster_pad_distance_{}_vs_{}'.format(cluster1+1, cluster2+1)
        #         h_title = 'Pad distance between: {} and {} clusters;N pads;N events'.format(cluster1+1, cluster2+1)
        #         self.h_dict[h_key] = TH1F(h_key, h_title, 400, 0, 60)

        #         # Distance in mm between cluster1 and cluster2
        #         h_key = 'h_cluster_coord_distance_{}_vs_{}'.format(cluster1+1, cluster2+1)
        #         h_title = 'Distance between: {} and {} clusters; d, [mm];N events'.format(cluster1+1, cluster2+1)
        #         self.h_dict[h_key] = TH1F(h_key, h_title, 500, 0, 100)

        #         # Energy ratio of cluster2 over cluster1
        #         h_key = 'h_cluster_ratio_{}_over_{}'.format(cluster2+1, cluster1+1)
        #         h_title = 'Energy ratio: {} over {} clusters;Ratio;N events'.format(cluster2+1, cluster1+1)
        #         self.h_dict[h_key] = TH1F(h_key, h_title, 400, 0, 1)

        #         # 2D: Distance in pads between clusters vs their energy ratio
        #         h_key = 'h_cluster_dist_ratio_for_{}_and_{}'.format(cluster1+1, cluster2+1)
        #         h_title = 'Distance vs ratio: {} and {} clusters;Distance;Ratio;N Events'.format(cluster1+1, cluster2+1)
        #         self.h_dict[h_key] = TH2F(h_key, '', 400, 0, 60, 400, 0, 1)

    def fill_histos(self):
        '''
        Fills all histograms (if created) with calculated
        variables or properties.
        '''
        # Number of clusters
        h_key = 'h_n_clusters'
        if h_key in self.h_dict:
            self.h_dict[h_key].Fill(len(self.cluster_list))

        # Energy of all 1pad clusters
        h_key = 'h_energy_1pad_clusters'
        if h_key in self.h_dict:
            for cluster in self.cluster_list:
                if cluster.n_pads == 1:
                    self.h_dict[h_key].Fill(cluster.energy)

        # Energy of all hits in tracker1
        h_key = 'h_tracker1_energy'
        if h_key in self.h_dict:
            self.h_dict[h_key].Fill(sum([signal.energy for signal in self.tracker1_data]))

        # Energy of all hits in tracker2
        h_key = 'h_tracker2_energy'
        if h_key in self.h_dict:
            self.h_dict[h_key].Fill(sum([signal.energy for signal in self.tracker2_data]))

        # 2D: Hits position weighted with energy in tracker1
        h_key = '2d_map_tracker1_energy'
        if h_key in self.h_dict:
            for signal in self.tracker1_data:
                self.h_dict[h_key].Fill(signal.sector, signal.pad, signal.energy)

        # 2D Hits position weighted with energy in tracker2
        h_key = '2d_map_tracker2_energy'
        if h_key in self.h_dict:
            for signal in self.tracker2_data:
                self.h_dict[h_key].Fill(signal.sector, signal.pad, signal.energy)

        # 2D: Hits positions in tracker1
        h_key = '2d_map_tracker1_hits'
        if h_key in self.h_dict:
            for signal in self.tracker1_data:
                self.h_dict[h_key].Fill(signal.sector, signal.pad)

        # 2D Hits positions in tracker2
        h_key = '2d_map_tracker2_hits'
        if h_key in self.h_dict:
            for signal in self.tracker2_data:
                self.h_dict[h_key].Fill(signal.sector, signal.pad)

        for idx, cluster in enumerate(self.cluster_list):
            if idx not in range(3):
                break
            h_key = 'h_cluster_energy_{}'.format(idx+1)
            if h_key in self.h_dict:
                self.h_dict[h_key].Fill(cluster.energy)

            h_key = 'h_cluster_pad_pos_{}'.format(idx+1)
            if h_key in self.h_dict:
                self.h_dict[h_key].Fill(cluster.pad)

            h_key = 'h_cluster_npads_{}'.format(idx+1)
            if h_key in self.h_dict:
                self.h_dict[h_key].Fill(cluster.n_pads)

        for idx1 in range(1):
            for idx2 in range(1, 3):
                if idx1 >= len(self.cluster_list) or idx2 >= len(self.cluster_list):
                    continue
                cluster1 = self.cluster_list[idx1]
                cluster2 = self.cluster_list[idx2]
                h_key = 'h_cluster_pad_distance_{}_vs_{}'.format(idx+1, idx2+1)
                if h_key in self.h_dict:
                    distance = abs(cluster1.pad-cluster2.pad)
                    self.h_dict[h_key].Fill(distance)

                h_key = 'h_cluster_coord_distance_{}_vs_{}'.format(idx+1, idx2+1)
                if h_key in self.h_dict:
                    distance = ((cluster1.x-cluster2.x)**2+(cluster1.y-cluster2.y)**2)**0.5
                    self.h_dict[h_key].Fill(distance)

                h_key = 'h_cluster_ratio_{}_over_{}'.format(idx2+1, idx+1)
                if h_key in self.h_dict:
                    energy_ratio = cluster2.energy/cluster1.energy
                    self.h_dict[h_key].Fill(energy_ratio)

                h_key = 'h_cluster_dist_ratio_for_{}_and_{}'.format(idx+1, idx2+1)
                if h_key in self.h_dict:
                    distance = abs(cluster1.pad-cluster2.pad)
                    self.h_dict[h_key].Fill(distance, energy_ratio)

        cluster1_pad = self.cluster_list[0].pad
        h_key = 'h_cluster_tracker1_distance'
        if h_key in self.h_dict:
            for signal in self.tracker1_data:
                distance = abs(signal.pad-cluster1_pad)
                self.h_dict[h_key].Fill(distance, signal.energy)

        h_key = 'h_cluster_tracker2_distance'
        if h_key in self.h_dict:
            for signal in self.tracker2_data:
                distance = abs(signal.pad-cluster1_pad)
                self.h_dict[h_key].Fill(distance, signal.energy)

        h_key = 'h_cluster_tracker1_position'
        if h_key in self.h_dict:
            for signal in self.tracker1_data:
                self.h_dict[h_key].Fill(signal.pad)

        h_key = 'h_cluster_tracker2_position'
        if h_key in self.h_dict:
            for signal in self.tracker2_data:
                self.h_dict[h_key].Fill(signal.pad)

    def extract_data(self, event):
        '''
        Extracts signals of event.THE MOST TIME CONSUMING FUNCTION :(
        '''
        # Make local variables
        position = self.position
        calib_energy = self.calib_energy

        # Read needed branches from the input ROOT file.
        id_arr = event.apv_id
        channel_arr = event.apv_ch
        signal_arr = event.apv_signal_maxfit
        apv_nn_output = event.apv_nn_output
        apv_fit_tau = event.apv_fit_tau
        apv_fit_t0 = event.apv_fit_t0
        apv_bint1 = event.apv_bint1
        # Lists for signals
        calorimeter_data = []
        tracker1_data = []
        tracker2_data = []

        # Loop through all signals(hits) in the event
        for hit in range(len(id_arr)):
            if (apv_fit_tau[hit] < 1 or apv_fit_tau[hit] > 3
               or signal_arr[hit] > 2000.
               or apv_fit_t0[hit] < (apv_bint1[hit]-2.7)
               or apv_fit_t0[hit] > (apv_bint1[hit]-0.5)):
                continue

            # Calculate position
            sector, pad, layer = position(id_arr[hit], channel_arr[hit])

            if (sector == 0 or sector == 3
               or (sector == 1 and pad < 20)
               or sector < 0  # This one is changed due to python C++ difference in %.
               or (layer < 2 and (signal_arr[hit] < 0. or apv_nn_output[hit] < 0.5))):
                continue

            energy = calib_energy(id_arr[hit], signal_arr[hit])
            # Ignore noisy cells in calorimeter
            if layer >= 2 and (energy < 1.4 or apv_nn_output[hit] < 0.5):
                continue

            # Choose what is data_list(tracker1/2,calorimeter)
            if layer == 0:
                data_list = tracker1_data
            elif layer == 1:
                data_list = tracker2_data
            else:
                data_list = calorimeter_data

            # If signal with this position already in the list: just add energy
            # Else: Add this signal to list and add energy
            for item in data_list:
                if (item.sector, item.pad) == (sector, pad):
                    item.energy += energy
                    break
            else:
                data_list.append(Signal(sector, pad, energy))

        # If no signals in calorimeter - skip event
        if len(calorimeter_data) == 0:
            return 0

        self.calorimeter_data = calorimeter_data
        self.tracker1_data = tracker1_data
        self.tracker2_data = tracker2_data
        # If everything is ok
        return 1

    def position(self, apv_id, apv_channel):
        '''
        Input: APV's id and channel
        Output:Tuple (sector, pad, layer): APV position in the detector
        Does: Read mapping array. Returns pad id and sector position of APV.
        Schematicaly pad ids and sectors numbering you can see in variables.py.
        '''

        # APV_id: odd - slave, even - master
        if apv_id < 4:
            if apv_id % 2 == 1:
                map_name = 'tb15_slave'
            else:
                map_name = 'tb15_master'
        elif apv_id >= 4 and apv_id < 14:
            if apv_id % 2 == 1:
                map_name = 'tb16_slave_divider'
            else:
                map_name = 'tb16_master_divider'
        elif apv_id == 14:
            map_name = 'tb16_master_tab_divider'
        elif apv_id == 15:
            map_name = 'tb16_slave_tab_divider'

        # Calculate corresponded position
        layer = apv_id//2
        # Case of -1 trace separetly. Because it caused a bug which is hard to find!!!
        # In C++ modulo "%" is different for negative numbers.
        # In python -1%64 result in 63. NOT -1 as expercted
        sector = apv_maps[map_name][apv_channel]//n_pads
        pad = apv_maps[map_name][apv_channel] % n_pads

        return sector, pad, layer

    def calib_energy(self, apv_id, apv_signal):
        '''
        Input:APV's id and signal(fit of RC-CR function) in MIPs.
        Output:Energy deposited in the layer
        Does: Reads calibration data files. Creates TGraphError with this data.
        Returns interpolated energy value - apv_energy
        '''
        # Optimization
        array = np.array

        # Calibration only to 1450 signal. No extrapolation. Just cut
        signal_treshold = 1450.
        cal_path = calib_path
        cal_file_name = calib_file_names[apv_id]

        # First point of callibraion curve
        x = [0.]
        y = [0.*16.5*1.164]
        x_err = [1.e-5]
        y_err = [1.e-5*16.5*1.164]

        # Calibration data in file written as (x,y,x_err,y_err) for each APV_id
        with open('%(cal_path)s%(cal_file_name)s' % locals(), 'r') as file:
            for i, line in enumerate(file):

                # skip a line with a title
                if i == 0:
                    continue

                # Calibration x-y data is inverted
                # Normalization D/MC according to L2 *1.09# /4.3 divide when No CD
                # *16.5*1.164 - is needed in order to get energy in MIPs.
                x.append(float(line.split('  ')[1]))
                y.append(float(line.split('  ')[0])*16.5*1.164)
                x_err.append(float(line.split('  ')[3]))
                y_err.append(float(line.split('  ')[2])*16.5*1.164)

        x = array(x)
        y = array(y)
        x_err = array(x_err)
        y_err = array(y_err)

        graph = TGraphErrors(len(x), x, y, x_err, y_err)

        # Take into account threshold
        if apv_signal > signal_treshold:
            signal = signal_treshold
        else:
            signal = apv_signal

        return graph.Eval(signal)

    def clustering_in_towers(self, merge='on'):
        '''
        Group signals into clusters changing their 'cluster' attribute.
        '''

        # Create list for cluster objects
        self.cluster_list = []

        # Find for every signal the most energetic neighbor
        # And write it as attribute 'neighbor'
        self.set_neighbors(self.calorimeter_data)

        # Collect signals into clusters, using linking-local-neighbor
        # algorithm.
        self.set_clusters(self.calorimeter_data)

        # Calculate number of primary clusters
        n_clusters = max([signal.cluster for signal in self.calorimeter_data])+1
        # Add cluster objects to the list.
        for cluster in range(n_clusters):
            self.cluster_list.append(Cluster(self.calorimeter_data, cluster))

        # If merge clusters option is on: merge clusters
        if merge == 'on':
            self.merge_clusters(self.calorimeter_data, self.cluster_list)

        # Analize only events with 1 cluster.
        if len(self.cluster_list) != 1:
            return 0
        # If everything is ok
        return 1

    def set_neighbors(self, data_list):
        '''
        Finds the most energetic neighbor among neighbors
        '''

        # Loop through all signals for which we are finding the neighbor:
        for signal in data_list:
            # define position of this signal
            center_sec, center_pad = signal.sector, signal.pad
            # create empty list for it's neighbors
            neighbors = []

            # loop over all neighbor signal-candidate
            for signal_neighbor in data_list:
                # If signal-candidate is in neighborhood to signal add it to neighbors list.
                # Neighbor can be signal itself! If it is the most energetic among it's neighbors
                if (signal_neighbor.sector in range(center_sec-1, center_sec+2)
                   and signal_neighbor.pad in range(center_pad-1, center_pad+2)):
                    neighbors.append(signal_neighbor)

            # Sort neighbors by energy
            neighbors_sorted = sorted(neighbors, key=lambda x: x.energy, reverse=True)

            # Pass the highest energetic neighbor to signal attribute "neighbor"
            signal.neighbor = neighbors_sorted[0]

    def set_clusters(self, data_list):
        '''
        Linking local neighbor algorithm.
        '''

        # Sort all signals by energy
        data_list = sorted(data_list, key=lambda x: x.energy, reverse=True)

        cluster_idx = 0
        # For signal in data list
        for signal in data_list:
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
            for signal in data_list:
                # If signal is not in a cluster: +1 to non_cluster counter
                # If the highest energetic neighbor of signal is in cluster
                # Add this signal to the same cluster
                if signal.cluster == -1:
                    n_non_clusters += 1
                    if signal.neighbor.cluster != -1:
                        signal.cluster = signal.neighbor.cluster

    def merge_clusters(self, data_list, cluster_list):
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
                        for item in data_list:
                            if item.cluster == cluster2_idx:
                                item.cluster = cluster1_idx
                        for item in data_list:
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


def langaus_fit(h_name):
    '''
    Fits energy distribution from ROOT file with Landau-Gaus
    convolution function.
    !!!!!ITS BROKEN. REPAIR !!!!!
    '''
    file = TFile('output.root', 'read')
    # Landau-Gauss fitting:
    fit_function = TF1('fit_function', langaufun, 0.1, 6, 4)
    fit_function.SetNpx(300)

    # Starting parameters
    fit_function.SetParameters(0.5, 1., 2600., 0.1)
    fit_function.SetParNames('Width', 'MP', 'Area', 'GSigma')

    print('It tries to fit. Please be pation and make yourself a tea :3')

    histo = file.Get(h_name)
    histo.SetTitle(h_name)

    histo.Fit('fit_function', "R")  # fit within specified range
    histo.Draw()

    print('Time used:', int(time.time()-start_time)//60, 'min ', end=' ')
    print(int(time.time()-start_time) % 60, 'sec')

    input('Pause. Enter a digit to exit')


def main():
    '''Analizing of data. Extraction/Clustering/Filling histos/Writing to file'''

    # Path/Name of input root file
    filename = './trees/run741_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root'
    # Create object of Class AnalizeCalorimeterEvent
    Analizer = EventData(filename)

    # Calcualte number of events to analize from Tree
    n_events = Analizer.tree.GetEntries()

    # Loop ove all events in tree
    for idx, event in enumerate(Analizer.tree):
        # For debuging. If you need to loop not over all events
        #if idx != 200:
        #    continue

        # Printout to see how many events are proceded/left.
        # And how much time spend per event.
        if idx % (500) == 0:
            time_min = (time.time()-start_time) // 60
            time_sec = (time.time()-start_time) % 60
            print('%(idx)i/%(n_events)i events' % locals(), end=' ')
            print('%(time_min)i min' % locals(), end=' ')
            print('%(time_sec)i sec' % locals())

        # Extract data.
        check = Analizer.extract_data(event)
        # If no hits - skip event
        if check == 0:
            continue

        # Do clustering to signals with merging
        check = Analizer.clustering_in_towers(merge='on')
        if check == 0:
            continue

        # Analize only events with 1 cluster
        filler = Analizer.fill_histos()

    # Update/create output file, where to write all histograms
    output_file = TFile('RENAME.root', 'update')

    # Write all histograms to the output root file
    for key in Analizer.h_dict.keys():
        Analizer.h_dict[key].Write()

    # Print Hooraay text
    input('Yaay I am finished :3')


# I believe It is faster if main code is written inside main().
main()
