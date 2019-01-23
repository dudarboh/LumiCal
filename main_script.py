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


class tower():
    '''
    Class to store information about signals in the event.
    It's position, energy, closest high-energy neighbor,
    cluster number of the signal
    '''

    def __init__(self, position, energy):
        self.position = position
        self.energy = energy
        self.neighbor = -1
        self.cluster = -1
        # Transform pad/sector numbers into x,y coordinates
        rho = 80+0.9+1.8*self.position[0]
        phi = np.pi/2+np.pi/12-np.pi/48-np.pi/24*position[1]
        self.coordinates = (rho*np.cos(phi), rho*np.sin(phi))


class AnalizeCalorimeterEvent(object):
    '''
    This class does calorimeter analysis:
    1) Extracts hits from imput root file (method extract_data)
    2.1) Gives every hit corresponding cluster number (method clustering_in_towers())
    2.2.) Merge some slucsters if some user's conditions are satisfied
    3) PlotCheck(), PlotClusterCheck() - plot 2d map of signals of 1 event.
    4) Fillxxx() - Fills histograms with corresponding values
    5) getxxx() - Compute corresponding value for the event
    6) position(), calib_energy() - calculate position and energy of the signal
    based on the input root file.
    7) merge_clusters() - merges all clusters.
    '''
    def __init__(self, filename):
        # Create root file and open "filename" file
        self.file = TFile.Open(filename)

        # Read a TTree from this file
        self.tree = self.file.apv_reco

        # Create dictionary to store all histograms
        self.h_dict = {}

    # Primary methods
    def extract_data(self, event):

        # List for signals data
        self.towers_list = []

        # Read all branches from the input ROOT file.
        id_arr = event.apv_id
        channel_arr = event.apv_ch
        signal_arr = event.apv_signal_maxfit
        apv_nn_output = event.apv_nn_output
        apv_fit_tau = event.apv_fit_tau
        apv_fit_t0 = event.apv_fit_t0
        apv_bint1 = event.apv_bint1

        # Make local variables to improve speed
        position = self.position
        calib_energy = self.calib_energy
        towers_list = self.towers_list

        # Loop through all signals(hits) in the event
        for hit in range(len(id_arr)):
            # Calculate sector, pad and layer(position) of the signal
            sector, pad, layer = position(id_arr[hit], channel_arr[hit])

            # Cuts. Analize: only calorimeter(layer>=2), only 2 central sectors, exclude bad mapping(pad<0)
            if (layer < 2
               or sector == 0 or sector == 3
               or (sector == 1 and pad < 20)
               or pad < 0
                # More cuts. it was ctrl+c ctrl+v from Sasha's code.
               or apv_fit_tau[hit] < 1 or apv_fit_tau[hit] > 3
               or signal_arr[hit] > 2000 or signal_arr[hit] < 0.75
               or apv_nn_output[hit] < 0.5
               or apv_fit_t0[hit] < (apv_bint1[hit]-2.7)
               or apv_fit_t0[hit] > (apv_bint1[hit]-0.5)):
                continue

            # Calculate energy of the event in MIPs
            energy = calib_energy(id_arr[hit], signal_arr[hit])

            # If bad event: skip event
            if energy < 0:
                return 0

            if energy < 1.4:
                continue

            # If signal for this sector,pad(NOT layer!) is in list: add energy
            # Else: Add this signal to list and add energy
            for item in towers_list:
                if item.position == (sector, pad):
                    item.energy += energy
                    break
            else:
                towers_list.append(tower((sector, pad), energy))

            # Write everything in global variable
            self.towers_list = towers_list

        # Return 1 if everything alright
        return 1

    def clustering_in_towers(self, merge='on'):

        # Make local variables to inprove speed.
        towers_list = self.towers_list
        merge_clusters = self.merge_clusters

        # for sll signals1 in the list:
        for signal in towers_list:
            # define position of this signal
            center_sec, center_pad = signal.position
            # create empty list for its neighbors
            neighbors = []

            # for all signals2 in the list
            for signal_neighbor in towers_list:
                # If signal2 is neighbor to signal1 (or signal2 == signal1) add it to neighbors list.
                # Signal2 can be signal1! If it is the most energetic among its neighbors
                if (signal_neighbor.position[0] in range(center_sec-1, center_sec+2)
                   and signal_neighbor.position[1] in range(center_pad-1, center_pad+2)):
                    neighbors.append(signal_neighbor)

            # Sort neighbors
            neighbors_sorted = sorted(neighbors, key=lambda x: x.energy, reverse=True)

            # Pass the highest energetic neighbor to signal1 attribute "neighbor"
            signal.neighbor = neighbors_sorted[0]

        # Sort all signals by energy
        towers_list = sorted(towers_list, key=lambda x: x.energy, reverse=True)

        cluster_idx = 0
        # For signal1 in data list
        for signal in towers_list:
            # If the neighbor of signal1 is signal1. (Means it local maximum)
            # Mark it as seed (give it cluster_index = 0, 1, 2, ...)
            if signal.neighbor.position == signal.position:
                signal.cluster = cluster_idx
                cluster_idx += 1

        # Variable to count number of non-cluster sigmals left
        n_non_clusters = -1
        # Stop when all signals got some cluster_index
        while n_non_clusters != 0:
            n_non_clusters = 0
            # for signal in data list
            for signal in towers_list:
                # If signal is not in a cluster: +1 to non_cluster counter
                # If the highest energetic neighbor of signal is in cluster
                # Add this signal to the same cluster
                if signal.cluster == -1:
                    n_non_clusters += 1
                    if signal.neighbor.cluster != -1:
                        signal.cluster = signal.neighbor.cluster

        # Write updated data to the global variable
        self.towers_list = towers_list

        # If merge clusters option is on:
        if merge == 'on':
            cluster1, cluster2 = 0, 0
            n_clusters = self.get_n_clusters()
            # Loops through all pairs of clusters
            while cluster1 in range(n_clusters):
                cluster2 = cluster1+1
                while cluster2 in range(n_clusters):
                    # For this pair: merge clusters()
                    merged = merge_clusters(cluster1, cluster2)
                    # If clusters were merged start looping from 0,0 again!
                    if merged:
                        cluster1, cluster2 = 0, 0
                    cluster2 += 1
                cluster1 += 1

    def merge_clusters(self, cluster1, cluster2):
        # Check, whether clusters cluster1 and cluster2 exist:
        ok1, ok2 = 0, 0
        for signal in self.towers_list:
            if signal.cluster == cluster1:
                ok1 = 1
            elif signal.cluster == cluster2:
                ok2 = 1
            if ok1 == 1 and ok2 == 1:
                break
        # If no: just skip merging
        else:
            return 0

        # Calculate distance and energy ration between clusters
        distance = self.get_pad_distance(cluster1, cluster2)
        ratio = self.get_energy_ratio(cluster1, cluster2)

        # Define if statement - when to merge clusters
        if (distance < 5
           or (distance < 20 and ratio < 0.1-0.1/20*distance)):
            # For all signals: if signal in cluster2: Make it cluster1.
            for signal in self.towers_list:
                if signal.cluster == cluster2:
                    signal.cluster = cluster1
                # Signals with cluster>cluster2: shift to the left (because cluster2 disappears).
                elif signal.cluster > cluster2:
                    signal.cluster -= 1
            # If were merged: return 1
            return 1
        else:
            return 0

    def PlotCheck(self, event):
        '''Plots map of signals in calorimeter towers as 2d histo.'''
        # Don't draw statistics on 2d map histo
        gStyle.SetOptStat(0)
        # Name of the histo
        h_key = 'check_event_{}'.format(event.apv_evt)
        # Local to improve speed (not necessary)
        h_dict = self.h_dict
        # Try is faster than if!
        # Fill histo. If does't exist: Create and fill
        try:
            for item in self.towers_list:
                h_dict[h_key].Fill(item.position[0], item.position[1], item.energy)
        except KeyError:
            h_dict[h_key] = TH2F(h_key, '', n_sectors+2, 0, n_sectors+2, n_pads, 0, n_pads)
            h_dict[h_key].SetTitle('event_{};sector;pad'.format(event.apv_evt))
            for item in self.towers_list:
                h_dict[h_key].Fill(item.position[0], item.position[1], item.energy)

        # Create Canvas and draw histo
        c1 = TCanvas('c1', h_key, 1800, 1800)
        h_dict[h_key].Draw("COLZTEXT")

        # Print picture to the png image
        c1.Print('./checks/'+h_key+'.png')
        # Plot statistic again. For others histos
        gStyle.SetOptStat(1)

    def PlotClusterCheck(self, event):
        '''Plots map of clusters as 2d histo'''
        # Inverse Palette, so 1st most energetic cluster will be the brightest
        TColor().InvertPalette()
        # Dont draw statistics
        gStyle.SetOptStat(0)

        # Name of the histo
        h_key = 'check_cluster_event_{}'.format(event.apv_evt)
        h_dict = self.h_dict
        # Try faster than if!
        # Try Fill, if KeyError: Create and Fill
        try:
            for item in self.towers_list:
                h_dict[h_key].Fill(item.position[0], item.position[1], item.cluster+1)
        except KeyError:
            h_dict[h_key] = TH2F(h_key, '', n_sectors+2, 0, n_sectors+2, n_pads, 0, n_pads)
            h_dict[h_key].SetTitle('cluster_event_{};sector;pad'.format(event.apv_evt))
            for item in self.towers_list:
                h_dict[h_key].Fill(item.position[0], item.position[1], item.cluster+1)

        # Create canvas and Draw
        c1 = TCanvas('c1', h_key, 1800, 1800)
        h_dict[h_key].Draw("COLZTEXT")

        # Print image and return palette and statistics settings on default values
        c1.Print('./checks/'+h_key+'.png')
        TColor().InvertPalette()
        gStyle.SetOptStat(1)

    def FillNclusters(self):
        '''Fill histo with number of clusters in event'''
        h_key = 'h_n_clusters'
        h_dict = self.h_dict
        try:
            h_dict[h_key].Fill(self.get_n_clusters())
        except KeyError:
            h_dict[h_key] = TH1F(h_key, '', 15, 0, 15)
            h_dict[h_key].SetTitle('N clusters;N Clusters;N Events')
            h_dict[h_key].Fill(self.get_n_clusters())

    def FillClusterEnergy(self, cluster):
        '''Fill histo with cluster energy in event'''
        # Check, whether this cluster number exist in this event
        for item in self.towers_list:
            if item.cluster == cluster:
                break
        else:
            return 0

        h_dict = self.h_dict
        h_key = 'h_energy_{}'.format(cluster+1)
        try:
            h_dict[h_key].Fill(self.get_cluster_energy(cluster))
        except KeyError:
            h_dict[h_key] = TH1F(h_key, '', 2000, 0, 500)
            h_dict[h_key].SetTitle('Energy: {} clust;Energy [MIP];N events'.format(cluster+1))
            h_dict[h_key].Fill(self.get_cluster_energy(cluster))

    def Fill1PadEnergy(self):
        '''Fills energy of all clusters with 1 pad size in the event'''
        h_dict = self.h_dict
        get_cluster_n_pads = self.get_cluster_n_pads
        get_cluster_energy = self.get_cluster_energy
        h_key = 'h_energy_1pad'
        try:
            for cluster in range(self.get_n_clusters()):
                if get_cluster_n_pads(cluster) == 1:
                    h_dict[h_key].Fill(get_cluster_energy(cluster))
        except KeyError:
            h_dict[h_key] = TH1F(h_key, '', 2000, 0, 20)
            h_dict[h_key].SetTitle('Energy: 1 pad clusters;Energy [MIP];N events')
            for cluster in range(self.get_n_clusters()):
                if get_cluster_n_pads(cluster) == 1:
                    h_dict[h_key].Fill(get_cluster_energy(cluster))

    def Fill1PadEnergy_dist_more_4_pads(self):
        '''
        Fills all clusters with size of 1 pad
        if they further than 4.5 pads from most energetic cluster
        '''
        h_dict = self.h_dict
        get_cluster_n_pads = self.get_cluster_n_pads
        get_cluster_energy = self.get_cluster_energy
        h_key = 'h_energy_1pad_dist_more_4_pads'
        try:
            for cluster in range(self.get_n_clusters()):
                if (get_cluster_n_pads(cluster) == 1
                   and self.get_pad_distance(1, cluster) > 4.5):
                    h_dict[h_key].Fill(get_cluster_energy(cluster))
        except KeyError:
            h_dict[h_key] = TH1F(h_key, '', 2000, 0, 20)
            h_dict[h_key].SetTitle('Energy: 1 pad clusters;Energy [MIP];N events')
            for cluster in range(self.get_n_clusters()):
                if (get_cluster_n_pads(cluster) == 1
                   and self.get_pad_distance(1, cluster) > 4.5):
                    h_dict[h_key].Fill(get_cluster_energy(cluster))

    def Fill1PadEnergy_dist_less_4_pads(self):
        '''
        Fills all clusters with size of 1 pad
        if they closer than 4.5 pads from most energetic cluster
        '''
        h_dict = self.h_dict
        get_cluster_n_pads = self.get_cluster_n_pads
        get_cluster_energy = self.get_cluster_energy
        h_key = 'h_energy_1pad_dist_less_4_pads'
        try:
            for cluster in range(self.get_n_clusters()):
                if (get_cluster_n_pads(cluster) == 1
                   and self.get_pad_distance(1, cluster) < 4.5):
                    h_dict[h_key].Fill(get_cluster_energy(cluster))
        except KeyError:
            h_dict[h_key] = TH1F(h_key, '', 2000, 0, 20)
            h_dict[h_key].SetTitle('Energy: 1 pad clusters;Energy [MIP];N events')
            for cluster in range(self.get_n_clusters()):
                if (get_cluster_n_pads(cluster) == 1
                   and self.get_pad_distance(1, cluster) < 4.5):
                    h_dict[h_key].Fill(get_cluster_energy(cluster))

    def FillClusterPadPos(self, cluster):
        '''Fill histo with pad position of the cluster in event'''
        for item in self.towers_list:
            if item.cluster == cluster:
                break
        else:
            return 0
        h_dict = self.h_dict
        h_key = 'h_position_{}'.format(cluster+1)
        try:
            h_dict[h_key].Fill(self.get_cluster_pad_pos(cluster))
        except KeyError:
            h_dict[h_key] = TH1F(h_key, '', 200, 0, n_pads)
            h_dict[h_key].SetTitle('Position: {} cluster;pos [pad];N events'.format(cluster+1))
            h_dict[h_key].Fill(self.get_cluster_pad_pos(cluster))

    def FillClusterNPads(self, cluster):
        '''Fills histo with number of pads of which consists cluster in event.'''
        for item in self.towers_list:
            if item.cluster == cluster:
                break
        else:
            return 0
        h_dict = self.h_dict
        h_key = 'h_npads_{}'.format(cluster+1)
        try:
            h_dict[h_key].Fill(self.get_cluster_n_pads(cluster))
        except KeyError:
            h_dict[h_key] = TH1F(h_key, '', 25, 0, 25)
            h_dict[h_key].SetTitle('N pads: {} cluster;N pads;N events'.format(cluster+1))
            h_dict[h_key].Fill(self.get_cluster_n_pads(cluster))

    def FillClusterPadDistance(self, cluster1, cluster2):
        '''Fill pad distance between cluster1 and cluster2 in event'''
        ok1, ok2 = 0, 0
        for item in self.towers_list:
            if item.cluster == cluster1:
                ok1 = 1
            elif item.cluster == cluster2:
                ok2 = 1
            if ok1 == 1 and ok2 == 1:
                break
        else:
            return 0
        h_dict = self.h_dict
        h_key = 'h_pad_distance_{}_vs_{}'.format(cluster1+1, cluster2+1)
        try:
            h_dict[h_key].Fill(self.get_pad_distance(cluster1, cluster2))
        except KeyError:
            h_dict[h_key] = TH1F(h_key, '', 400, 0, 60)
            self.h_dict[h_key].SetTitle('Pad distance between: {} and {} clusters;N pads;\
                                            N events'.format(cluster1+1, cluster2+1))
            h_dict[h_key].Fill(self.get_pad_distance(cluster1, cluster2))

    def FillClusterCoordDistance(self, cluster1, cluster2):
        '''Fill pad distance between cluster1 and cluster2 in event'''
        ok1, ok2 = 0, 0
        for item in self.towers_list:
            if item.cluster == cluster1:
                ok1 = 1
            elif item.cluster == cluster2:
                ok2 = 1
            if ok1 == 1 and ok2 == 1:
                break
        else:
            return 0
        h_dict = self.h_dict
        h_key = 'h_coord_distance_{}_vs_{}'.format(cluster1+1, cluster2+1)
        try:
            h_dict[h_key].Fill(self.get_coord_distance(cluster1, cluster2))
        except KeyError:
            h_dict[h_key] = TH1F(h_key, '', 400, 0, 60)
            self.h_dict[h_key].SetTitle('Distance between: {} and {} clusters; d, [mm];\
                                            N events'.format(cluster1+1, cluster2+1))
            h_dict[h_key].Fill(self.get_coord_distance(cluster1, cluster2))

    def FillClusterRatio(self, cluster1, cluster2):
        '''Fill histo with ratio of cluster energies in event'''
        ok1, ok2 = 0, 0
        for item in self.towers_list:
            if item.cluster == cluster1:
                ok1 = 1
            elif item.cluster == cluster2:
                ok2 = 1
            if ok1 == 1 and ok2 == 1:
                break
        else:
            return 0
        h_dict = self.h_dict
        h_key = 'h_ratio_{}_over_{}'.format(cluster2+1, cluster1+1)
        try:
            h_dict[h_key].Fill(self.get_energy_ratio(cluster1, cluster2))
        except KeyError:
            h_dict[h_key] = TH1F(h_key, '', 400, 0, 1)
            h_dict[h_key].SetTitle('Energy ratio: {} over {} clusters;Ratio;\
                                    N events'.format(cluster2+1, cluster1+1))
            h_dict[h_key].Fill(self.get_energy_ratio(cluster1, cluster2))

    def FillClusterInverseRatio(self, cluster1, cluster2):
        '''Fill ratio of clusters energy, but higher_energy/lower_energy (inversed)'''
        ok1, ok2 = 0, 0
        for item in self.towers_list:
            if item.cluster == cluster1:
                ok1 = 1
            elif item.cluster == cluster2:
                ok2 = 1
            if ok1 == 1 and ok2 == 1:
                break
        else:
            return 0
        h_dict = self.h_dict
        h_key = 'h__inverse_ratio_{}_over_{}'.format(cluster2+1, cluster1+1)
        try:
            h_dict[h_key].Fill(1/self.get_energy_ratio(cluster1, cluster2))
        except KeyError:
            h_dict[h_key] = TH1F(h_key, '', 400, 0, 200)
            h_dict[h_key].SetTitle('Inverse Energy ratio: {} over {} clusters;Ratio;\
                                    N events'.format(cluster2+1, cluster1+1))
            h_dict[h_key].Fill(1/self.get_energy_ratio(cluster1, cluster2))

    def FillClusterDistVsRatio(self, cluster1, cluster2):
        '''Fills 2d histogram of energy ratio vs distance for cluster1 vs cluster2 in event'''
        ok1, ok2 = 0, 0
        for item in self.towers_list:
            if item.cluster == cluster1:
                ok1 = 1
            elif item.cluster == cluster2:
                ok2 = 1
            if ok1 == 1 and ok2 == 1:
                break
        else:
            return 0
        h_dict = self.h_dict
        h_key = 'h_dist_ratio_for_{}_and_{}'.format(cluster1+1, cluster2+1)
        try:
            h_dict[h_key].Fill(self.get_pad_distance(cluster1, cluster2),
                               self.get_energy_ratio(cluster1, cluster2))
        except KeyError:
            h_dict[h_key] = TH2F(h_key, '', 400, 0, 60, 400, 0, 1)
            h_dict[h_key].SetTitle('Distance vs ratio: {} and {} clusters;Distance;\
                                    Ratio;N Events'.format(cluster1+1, cluster2+1))
            h_dict[h_key].Fill(self.get_pad_distance(cluster1, cluster2),
                               self.get_energy_ratio(cluster1, cluster2))

    # Secondary methods
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
        sector = apv_maps[map_name][apv_channel]//n_pads
        pad = apv_maps[map_name][apv_channel] % n_pads
        layer = apv_id//2

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
                x.append(float(line.split('  ')[1]))
                y.append(float(line.split('  ')[0])*16.5*1.164)
                x_err.append(float(line.split('  ')[3]))
                y_err.append(float(line.split('  ')[2])*16.5*1.164)

        x = array(x)
        y = array(y)
        x_err = array(x_err)
        y_err = array(y_err)

        graph = TGraphErrors(len(x), x, y, x_err, y_err)

        # Copypasted from Sasha's code. Scale calibration graphs
        # Normalization D/MC according to L2 *1.09# /4.3 divide when No CD
        # for point in range(graph.GetN()):
            # Scale Y to get MIPS
        #    graph.GetY()[point] *= 16.5*1.164

            # y_err zero anyway
        #    graph.GetEY()[point] *= 16.5*1.164

        # Take into account threshold
        if apv_signal > signal_treshold:
            signal = signal_treshold
        else:
            signal = apv_signal

        return graph.Eval(signal)

    # Get functions

    def get_n_clusters(self):
        '''Calculate number of clusters in event'''
        # List of cluster_indices of data list
        cluster_list = [item.cluster for item in self.towers_list]
        # If no data points, return 0
        if not cluster_list:
            return 0
        # Return maximum cluster_index+1
        return max(cluster_list)+1

    def get_cluster_energy(self, cluster):
        '''Calculates cluster energy'''
        # List of energy of all data points which have attribute cluster == cluster
        energy_list = [item.energy for item in self.towers_list if item.cluster == cluster]
        # Return sum of all these energies
        return sum(energy_list)

    def get_cluster_pad_pos(self, cluster):
        '''Calculates pad position of the cluster'''
        cluster_pos = 0
        # Calculate energy of the cluster
        cluster_energy = self.get_cluster_energy(cluster)
        # Create list of tuples(pad position, energy) for all signals of cluster
        pos_energy_list = [(item.position[1], item.energy) for item in self.towers_list if item.cluster == cluster]

        # Calculate position as sum with weights(energies) for each point
        for pos_energy in pos_energy_list:
            cluster_pos += pos_energy[0]*pos_energy[1]/cluster_energy
        # if no cluster returns 0 pos
        return cluster_pos

    def get_cluster_sector_pos(self, cluster):
        '''Calculates sector position of the cluster'''
        cluster_pos = 0
        # Calculate cluster energy
        cluster_energy = self.get_cluster_energy(cluster)
        # Create a list of tuples(sector position, energy) of each signal of cluster
        pos_energy_list = [(item.position[0], item.energy) for item in self.towers_list if item.cluster == cluster]

        # Calculate position as sum with weights(energies) over all points
        for pos_energy in pos_energy_list:
            cluster_pos += pos_energy[0]*pos_energy[1]/cluster_energy
        # if no cluster returns 0 pos
        return cluster_pos

    def get_cluster_coord(self, cluster):
        '''Calculates x,y position of the cluster'''
        cluster_x = 0
        cluster_y = 0
        # Calculate energy of the cluster
        cluster_energy = self.get_cluster_energy(cluster)
        # Create list of tuples(pad position, energy) for all signals of cluster
        x_energy_list = [(item.coordinates[0], item.energy) for item in self.towers_list if item.cluster == cluster]

        y_energy_list = [(item.coordinates[1], item.energy) for item in self.towers_list if item.cluster == cluster]

        # Calculate position as sum with weights(energies) for each point
        for x in x_energy_list:
            cluster_x += x[0]*x[1]/cluster_energy
        for y in y_energy_list:
            cluster_y += y[0]*y[1]/cluster_energy

        # if no cluster returns 0 coordinates
        return cluster_x, cluster_y

    def get_cluster_n_pads(self, cluster):
        '''Calculates number of pads(towers) in cluster'''
        # Create list of all signals(towers) in cluster
        n_pads_list = [item for item in self.towers_list if item.cluster == cluster]

        # Return length of the list(number of towers).
        return len(n_pads_list)

    def get_pad_distance(self, cluster1, cluster2):
        '''Calculate pad(projected) distance between 2 clusters'''

        # Just for speed performance (not necessary)
        get_cluster_pad_pos = self.get_cluster_pad_pos

        # Calculate position of the clusters and return its difference
        position1 = get_cluster_pad_pos(cluster1)
        position2 = get_cluster_pad_pos(cluster2)
        # This is only projection on pad distance!
        return abs(position1-position2)

    def get_sector_distance(self, cluster1, cluster2):
        '''Calculate sector(projected) distance between 2 clusters'''

        # Just for speed (not necessary)
        get_cluster_sector_pos = self.get_cluster_sector_pos

        # Calculate sector positions and return difference
        position1 = get_cluster_sector_pos(cluster1)
        position2 = get_cluster_sector_pos(cluster2)
        # This is only projection on sector distance!
        return abs(position1-position2)

    def get_coord_distance(self, cluster1, cluster2):
        '''Calculate x,y distance between 2 clusters'''

        # Calculate sector positions and return difference
        x1, y1 = self.get_cluster_coord(cluster1)
        x2, y2 = self.get_cluster_coord(cluster2)
        # This is only projection on sector distance!
        return np.sqrt((x1-x2)**2+(y1-y2)**2)

    def get_energy_ratio(self, cluster1, cluster2):
        '''Calculates ratio of energy between 2 clusters'''

        # Just for performance
        get_cluster_energy = self.get_cluster_energy

        # Calculates enegies of the clusters and return their ratio
        energy1 = get_cluster_energy(cluster1)
        energy2 = get_cluster_energy(cluster2)
        return energy2/energy1


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
    Analizer = AnalizeCalorimeterEvent(filename)

    # Calcualte number of events to analize from Tree
    n_events = Analizer.tree.GetEntries()

    # Loop ove all events in tree
    for idx, event in enumerate(Analizer.tree):
        # For debuging. If you need to loop not over all events
        # if idx == 1000:
        #     break

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
        # if bad data(energy<0) - skip event
        if check == 0:
            continue

        # Do clustering to signals with merging
        Analizer.clustering_in_towers(merge='on')

        # Fill according histograms. For details open FillName() methods
        Analizer.FillNclusters()
        for cluster in range(1, 4):
            Analizer.FillClusterCoordDistance(0, cluster)
            Analizer.FillClusterRatio(0, cluster)
            Analizer.FillClusterDistVsRatio(0, cluster)

        for cluster in range(0, 4):
            Analizer.FillClusterEnergy(cluster)
            Analizer.FillClusterPadPos(cluster)
            Analizer.FillClusterNPads(cluster)

    # Update/create output file, where to write all histograms
    output_file = TFile('RENAME.root', 'update')

    # Write all histograms to the output root file
    for key in Analizer.h_dict.keys():
        Analizer.h_dict[key].Write()

    # Print Hooraay text
    input('Yaay I am finished :3')


# I believe It is faster if main code is written inside main().
main()
