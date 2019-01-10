'''
This file does analysis of cluster reconstruction in calorimeter and trackers
signals in coresponded functions cluster_analysis(), tracker_analysis().
To start this script user needs all imported libraries. Python3 with root.
P.S. It may work with python2. But it may HAVE BUGS!!! Because in python2:
1/2=0, when in python3 1/2=0.5 and so on... + python2 slower x4 times. So adapt
to use python3. You will need ROOT with python3. If it works only with python2
but not python3. You have to reinstall it and compile source (dont use binary)
with python3 executable/library/include paths provided in ./configure stage.
Include directory for me was tricky: /path/to/python/include/python3.7m. So
play
round it.
'''

from ROOT import TH1F, TH2F, TF1, TGraphErrors, TFile, gROOT, TCanvas, gStyle
import time
import numpy as np

from variables import n_sectors, n_pads, apv_maps, calib_path, \
    calib_file_names, langaufun

start_time = time.time()

gROOT.SetBatch(1)
gStyle.SetOptStat(0)


class AnalizeCalorimeterEvent(object):
    def __init__(self, filename):
        self.file = TFile.Open(filename)
        self.tree = self.file.apv_reco
        self.h_dict = {}

    # Special methods
    def __getitem__(self, key):
        return getattr(self, key)

    # Primary methods
    def extract_data(self, event):

        # Create data list
        # Full data
        # self.data_arr = np.zeros((4, 64, 8))
        # Towers data
        self.towers_arr = np.zeros((4, 64))

        id_arr = event.apv_id
        channel_arr = event.apv_ch
        signal_arr = event.apv_signal_maxfit

        for hit in range(len(id_arr)):
            sector, pad, layer = self.position(id_arr[hit], channel_arr[hit])

            # Analize only calorimeter
            if layer < 2:
                continue
            # Geometry cuts
            if not (sector == 1 or sector == 2):
                continue
            # Cut on APV noisy area
            if sector == 1 and pad < 20:
                continue
            # Cut on APV maping
            if pad < 0:
                continue
            # More cuts. it was ctrl+c ctrl+v from Sasha's code.
            cond_1 = event.apv_fit_tau[hit] > 1 and event.apv_fit_tau[hit] < 3
            cond_2 = signal_arr[hit] < 2000
            cond_31 = layer < 2 and signal_arr[hit] > 0
            cond_32 = event.apv_nn_output[hit] > 0.5
            cond_41 = layer >= 2 and signal_arr[hit] > 0.75
            cond_42 = event.apv_nn_output[hit] > 0.5
            cond_51 = event.apv_fit_t0[hit] > (event.apv_bint1[hit]-2.7)
            cond_52 = event.apv_fit_t0[hit] < (event.apv_bint1[hit]-0.5)

            if not (cond_1 and cond_2
                    and ((cond_31 and cond_32) or (cond_41 and cond_42))
                    and cond_51 and cond_52):
                continue

            energy = self.calib_energy(id_arr[hit], signal_arr[hit])
            # self.data_arr[sector, pad, layer] += energy
            self.towers_arr[sector, pad] += energy

    def clustering_in_towers(self, merging='on'):

        neighbors_arr = np.full((4, 64, 2), -1)
        self.clusters_arr = np.full((4, 64), -1)

        for sec in range(n_sectors):
            for pad in range(n_pads):
                max_energy = 0
                for s_itr in [-1, 0, 1]:
                    for p_itr in [-1, 0, 1]:
                        if not (0 <= (sec+s_itr) < n_sectors
                                and 0 <= (pad+p_itr) < n_pads
                                and self.towers_arr[sec, pad] != 0):
                            continue
                        if max_energy > self.towers_arr[sec+s_itr, pad+p_itr]:
                            continue
                        neighbors_arr[sec, pad] = (sec+s_itr, pad+p_itr)
                        max_energy = self.towers_arr[sec+s_itr, pad+p_itr]

        cluster_idx = 0
        start_energy = np.amax(self.towers_arr[self.clusters_arr == -1])

        while start_energy > 0:
            hit_idx = np.argwhere(self.towers_arr == start_energy)[0]
            self.collect_cluster(hit_idx, neighbors_arr, cluster_idx)
            cluster_idx += 1
            start_energy = np.amax(self.towers_arr[self.clusters_arr == -1])
            if cluster_idx > 1000:
                print('WEEEEEEEIRD EVENT!!!!!!!!!!!29441???')
                break
        if merging == 'on':
            for cluster in range(np.amax(self.clusters_arr)):
                self.merge_clusters(cluster, cluster+1)

    def merge_clusters(self, cluster1, cluster2):
        if not ((self.clusters_arr == cluster1).any()
                and (self.clusters_arr == cluster2).any()):
            return 0
        elif (self.get_sector_distance(cluster1, cluster2) < 1.5
              and self.get_pad_distance(cluster1, cluster2) < 4.5):
                # Make 0 cluster 1st, and then substarct 1
                self.clusters_arr[self.clusters_arr == cluster1] = cluster2
                for sec in range(n_sectors):
                    for pad in range(n_pads):
                        if self.clusters_arr[sec, pad] >= cluster2:
                            self.clusters_arr[sec, pad] -= 1
                self.merge_clusters(cluster1, cluster2)
        else:
            self.merge_clusters(cluster1, cluster2+1)

    def PlotCheck(self, event):
        h_key = 'check_event_{}'.format(event.apv_evt)
        if h_key not in self.h_dict.keys():
            self.h_dict[h_key] = TH2F(h_key, '', 6, 0, 6, 64, 0, 64)
            self.h_dict[h_key].SetTitle('event_{};sector;pad'.format(event.apv_evt))
        for sec in range(n_sectors):
            for pad in range(n_pads):
                self.h_dict[h_key].Fill(sec, pad, self.towers_arr[sec, pad])
        c1 = TCanvas('c1', h_key, 1800, 1800)
        self.h_dict[h_key].Draw("COLZ")

        c1.Print('./check_pics/'+h_key+'.png')

    def FillNclusters(self):
        h_key = 'h_n_clusters'
        if h_key not in self.h_dict.keys():
            self.h_dict[h_key] = TH1F(h_key, '', 15, 0, 15)
            self.h_dict[h_key].SetTitle('N clusters;N Clusters;N Events')
        if (self.towers_arr >= 0).all():
            self.h_dict[h_key].Fill(self.get_n_clusters())

    def FillClusterEnergy(self, cluster):
        if not (self.clusters_arr == cluster).any():
            return 0
        h_key = 'h_energy_{}'.format(cluster+1)
        if h_key not in self.h_dict.keys():
            self.h_dict[h_key] = TH1F(h_key, '', 2000, 0, 500)
            title = 'Energy: {} clust;Energy [MIP];N events'.format(cluster+1)
            self.h_dict[h_key].SetTitle(title)
        if (self.towers_arr >= 0).all():
            self.h_dict[h_key].Fill(self.get_cluster_energy(cluster))

    def Fill1PadEnergy(self):
        h_key = 'h_energy_1pad'
        if h_key not in self.h_dict.keys():
            self.h_dict[h_key] = TH1F(h_key, '', 2000, 0, 500)
            title = 'Energy: 1 pad clusters;Energy [MIP];N events'
            self.h_dict[h_key].SetTitle(title)
        if (self.towers_arr >= 0).all():
            for cluster in range(self.get_n_clusters()):
                if self.get_cluster_n_pads(cluster) == 1:
                    self.h_dict[h_key].Fill(self.get_cluster_energy(cluster))

    def FillClusterPadPos(self, cluster):
        if not (self.clusters_arr == cluster).any():
            return 0
        h_key = 'h_position_{}'.format(cluster+1)
        if h_key not in self.h_dict.keys():
            self.h_dict[h_key] = TH1F(h_key, '', 200, 0, 64)
            title = 'Position: {} cluster;pos [pad];N events'.format(cluster+1)
            self.h_dict[h_key].SetTitle(title)
        if (self.towers_arr >= 0).all():
            self.h_dict[h_key].Fill(self.get_cluster_pad_pos(cluster))

    def FillClusterNPads(self, cluster):
        if not (self.clusters_arr == cluster).any():
            return 0
        h_key = 'h_npads_{}'.format(cluster+1)
        if h_key not in self.h_dict.keys():
            self.h_dict[h_key] = TH1F(h_key, '', 25, 0, 25)
            title = 'N pads: {} cluster;N pads;N events'.format(cluster+1)
            self.h_dict[h_key].SetTitle(title)
        if (self.towers_arr >= 0).all():
            self.h_dict[h_key].Fill(self.get_cluster_n_pads(cluster))

    def FillClusterDistance(self, cluster1, cluster2):
        if not ((self.clusters_arr == cluster1).any()
                and (self.clusters_arr == cluster2).any()):
            return 0
        h_key = 'h_distance_{}_vs_{}'.format(cluster1+1, cluster2+1)
        if h_key not in self.h_dict.keys():
            self.h_dict[h_key] = TH1F(h_key, '', 400, 0, 60)
            title = 'Distance between: {} and {} clusters;N pads;\
                     N events'.format(cluster1+1, cluster2+1)
            self.h_dict[h_key].SetTitle(title)
        if (self.towers_arr >= 0).all():
            self.h_dict[h_key].Fill(self.get_pad_distance(cluster1, cluster2))

    def FillClusterRatio(self, cluster1, cluster2):
        if not ((self.clusters_arr == cluster1).any()
                and (self.clusters_arr == cluster2).any()):
            return 0
        h_key = 'h_ratio_{}_over_{}'.format(cluster2+1, cluster1+1)
        if h_key not in self.h_dict.keys():
            self.h_dict[h_key] = TH1F(h_key, '', 400, 0, 1)
            title = 'Energy ratio: {} over {} clusters;Ratio;\
                     N events'.format(cluster2+1, cluster1+1)
            self.h_dict[h_key].SetTitle(title)
        if (self.towers_arr >= 0).all():
            self.h_dict[h_key].Fill(self.get_energy_ratio(cluster1, cluster2))

    def FillClusterInverseRatio(self, cluster1, cluster2):
        if not ((self.clusters_arr == cluster1).any()
                and (self.clusters_arr == cluster2).any()):
            return 0
        h_key = 'h__inverse_ratio_{}_over_{}'.format(cluster2+1, cluster1+1)
        if h_key not in self.h_dict.keys():
            self.h_dict[h_key] = TH1F(h_key, '', 400, 0, 200)
            title = 'Inverse Energy ratio: {} over {} clusters;Ratio;\
                     N events'.format(cluster2+1, cluster1+1)
            self.h_dict[h_key].SetTitle(title)
        if (self.towers_arr >= 0).all():
            self.h_dict[h_key].Fill(1/self.get_energy_ratio(cluster1, cluster2))

    def FillClusterDistVsRatio(self, cluster1, cluster2):
        if not ((self.clusters_arr == cluster1).any()
                and (self.clusters_arr == cluster2).any()):
            return 0
        h_key = 'h_dist_ratio_for_{}_and_{}'.format(cluster1+1, cluster2+1)
        if h_key not in self.h_dict.keys():
            self.h_dict[h_key] = TH2F(h_key, '', 400, 0, 60, 400, 0, 1)
            title = 'Distance vs ratio: {} and {} clusters;Distance;\
                     Ratio;N Events'.format(cluster1+1, cluster2+1)
            self.h_dict[h_key].SetTitle(title)
        if (self.towers_arr >= 0).all():
            self.h_dict[h_key].Fill(self.get_pad_distance(cluster1, cluster2),
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

        # Calibration only to 1450 signal. No extrapolation. Just cut
        signal_treshold = 1450.

        # First point of callibraion curve
        x = [0.]
        y = [0.]
        x_err = [1.e-5]
        y_err = [1.e-5]

        # Calibration data in file written as (x,y,x_err,y_err) for each APV_id
        with open(calib_path+calib_file_names[apv_id], 'r') as file:
            for i, line in enumerate(file):

                # skip a line with a title
                if i == 0:
                    continue

                # Calibration x-y data is inverted
                x += [float(line.split('  ')[1])]
                y += [float(line.split('  ')[0])]
                x_err += [float(line.split('  ')[3])]
                y_err += [float(line.split('  ')[2])]

        x = np.array(x)
        y = np.array(y)
        x_err = np.array(x_err)
        y_err = np.array(y_err)

        graph = TGraphErrors(len(x), x, y, x_err, y_err)

        # Copypasted from Sasha's code. Scale calibration graphs
        # Normalization D/MC according to L2 *1.09# /4.3 divide when No CD
        for point in range(graph.GetN()):
            # Scale Y to get MIPS
            graph.GetY()[point] *= 16.5*1.164

            # y_err zero anyway
            graph.GetEY()[point] *= 16.5*1.164

        # Take into account threshold
        if apv_signal > signal_treshold:
            signal = signal_treshold
        else:
            signal = apv_signal

        return graph.Eval(signal)

    def collect_cluster(self, hit_idx, neighbors, cluster_idx):
        self.clusters_arr[hit_idx[0], hit_idx[1]] = cluster_idx
        for sec in range(n_sectors):
            for pad in range(n_pads):
                if (np.all(neighbors[sec, pad] == hit_idx)
                   and self.clusters_arr[sec, pad] == -1):
                    self.collect_cluster((sec, pad), neighbors, cluster_idx)

    # Get functions

    def get_n_clusters(self):
        if (self.towers_arr >= 0).all():
            return np.amax(self.clusters_arr)+1
        else:
            return -1000

    def get_cluster_energy(self, cluster):
        if not ((self.clusters_arr == cluster).any()
                and (self.towers_arr >= 0).all()):
            return -1000
        return np.sum(self.towers_arr[self.clusters_arr == cluster])

    def get_cluster_pad_pos(self, cluster):
        if not ((self.clusters_arr == cluster).any()
                and (self.towers_arr >= 0).all()):
            return -1000
        cluster_pos = 0
        cluster_energy = self.get_cluster_energy(cluster)
        for sec in range(n_sectors):
            for pad in range(n_pads):
                if self.clusters_arr[sec, pad] != cluster:
                    continue
                cluster_pos += pad*self.towers_arr[sec, pad]/cluster_energy
        return cluster_pos

    def get_cluster_sector_pos(self, cluster):
        if not ((self.clusters_arr == cluster).any()
                and (self.towers_arr >= 0).all()):
            return -1000
        cluster_pos = 0
        cluster_energy = self.get_cluster_energy(cluster)
        for sec in range(n_sectors):
            for pad in range(n_pads):
                if self.clusters_arr[sec, pad] != cluster:
                    continue
                cluster_pos += sec*self.towers_arr[sec, pad]/cluster_energy
        return cluster_pos

    def get_cluster_n_pads(self, cluster):
        if not ((self.clusters_arr == cluster).any()
                and (self.towers_arr >= 0).all()):
            return -1000
        return np.sum(self.clusters_arr == cluster)

    def get_pad_distance(self, cluster1, cluster2):
        if not ((self.clusters_arr == cluster1).any()
                and (self.clusters_arr == cluster2).any()
                and (self.towers_arr >= 0).all()):
            return -1000
        position1 = self.get_cluster_pad_pos(cluster1)
        position2 = self.get_cluster_pad_pos(cluster2)
        # This is only projection on pad distance!
        return abs(position1-position2)

    def get_sector_distance(self, cluster1, cluster2):
        if not ((self.clusters_arr == cluster1).any()
                and (self.clusters_arr == cluster2).any()
                and (self.towers_arr >= 0).all()):
            return -1000
        position1 = self.get_cluster_sector_pos(cluster1)
        position2 = self.get_cluster_sector_pos(cluster2)
        # This is only projection on pad distance!
        return abs(position1-position2)

    def get_energy_ratio(self, cluster1, cluster2):
        if not ((self.clusters_arr == cluster1).any()
                and (self.clusters_arr == cluster2).any()
                and (self.towers_arr >= 0).all()):
            return -1000
        energy1 = self.get_cluster_energy(cluster1)
        energy2 = self.get_cluster_energy(cluster2)
        return energy2/energy1


def tracker_analysis(tracker_layer):
    '''
    Calculates next variables for a tracker layer: Number of hits,
    total energy deposited, positions of hits. Plots histograms.
    '''

    gROOT.SetBatch(1)
    # gStyle.SetOptStat(0)

    input_filename = 'run741_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root'
    input_file = TFile.Open(input_filename)

    h_energy_calib = TH1F('h_energy_calib', '', 100, 0, 10.)
    h_energy_calib.SetTitle('Total energy per event;Energy, [MIP];N events')

    h_position = TH1F('h_position', '', 64, 0-0.5, 64-0.5)
    h_position.SetTitle('Position per hit;Pad number;N events')

    h_n_hits = TH1F('h_hits_number', '', 5, 0-0.5, 5-0.5)
    h_n_hits.SetTitle('N Hits per event;N hits;N events')

    for idx, event in enumerate(input_file.apv_reco):
        # Use to debug to pass not all 50k points
        if idx == 200:
            break

        n_events = input_file.apv_reco.GetEntries()
        if idx == 0:
            print('Evaluating', n_events, 'Events:')

        if idx % (n_events*2//100) == 0:
            print(tracker_layer, 'Tracker:', end=' ')
            print(100*idx//n_events, '% are done', end=' ')
            print('Time used:', (time.time()-start_time)//60, 'min', end=' ')
            print(int((time.time()-start_time)) % 60, 'sec')

        n_hits = 0
        energy = 0.

        id_arr = event.apv_id
        channel_arr = event.apv_ch
        signal_arr = event.apv_signal_maxfit

        for j in range(len(id_arr)):

            sector, pad, layer = position(id_arr[j], channel_arr[j])

            # ###Start of cut section###
            # Cut on APV maping
            if pad < 0:
                continue
            # Geometry cuts
            if not (sector == 1 or sector == 2):
                continue
            # Cut on APV noisy area
            if sector == 1 and pad < 20:
                continue

            # Analyse only defined trackers
            if layer != tracker_layer:
                continue

            # More cuts. it was ctrl+c ctrl+v from Sasha's code.
            cond_1 = event.apv_fit_tau[j] > 1. and event.apv_fit_tau[j] < 3.
            cond_2 = signal_arr[j] < 2000.
            cond_31 = layer < 2 and signal_arr[j] > 0.
            cond_32 = event.apv_nn_output[j] > 0.5
            cond_41 = layer >= 2 and signal_arr[j] > 0.75
            cond_42 = event.apv_nn_output[j] > 0.5
            cond_51 = event.apv_fit_t0[j] > (event.apv_bint1[j]-2.7)
            cond_52 = event.apv_fit_t0[j] < (event.apv_bint1[j]-0.5)

            if not (cond_1 and cond_2
               and ((cond_31 and cond_32) or (cond_41 and cond_42))
               and cond_51 and cond_52):
                continue

            # ###End of cut section###

            n_hits += 1
            energy += calib_energy(id_arr[j], signal_arr[j])

            h_position.Fill(pad)

        h_energy_calib.Fill(energy)
        h_n_hits.Fill(n_hits)

    output_file = TFile('output.root', 'update')

    h_energy_calib.Write('energy_tracker_'+str(tracker_layer))
    h_position.Write('position_tracker_'+str(tracker_layer))
    h_n_hits.Write('hits_tracker_'+str(tracker_layer))

    output_file.Close()


def energy_fit(tracker_layer):
    file = TFile('output.root', 'read')
    # Landau-Gauss fitting:
    fit_function = TF1('fit_function', langaufun, 0.1, 6, 4)
    fit_function.SetNpx(300)

    # Starting parameters
    fit_function.SetParameters(0.5, 1., 2600., 0.1)
    fit_function.SetParNames('Width', 'MP', 'Area', 'GSigma')

    print('It tries to fit. Please be pation and make yourself a tea :3')

    histo = file.Get('energy_tracker_'+str(tracker_layer))
    histo.SetTitle('Total Energy per event: Tracker '+str(tracker_layer))

    histo.Fit('fit_function', "R")  # fit within specified range
    histo.Draw()

    print('Time used:', int(time.time()-start_time)//60, 'min ', end=' ')
    print(int(time.time()-start_time) % 60, 'sec')

    input('Pause. Enter a digit to exit')


filename = './trees/run741_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root'
Analizer = AnalizeCalorimeterEvent(filename)
n_events = Analizer.tree.GetEntries()

for idx, event in enumerate(Analizer.tree):
    #if idx == 200:
    #    break

    if idx == 0:
        print(n_events, 'events')
    elif idx % (n_events*0.1//100) == 0:
        print(idx, 'event', end=' ')
        print(int(time.time()-start_time)//60, 'min', end=' ')
        print(int((time.time()-start_time)) % 60, 'sec')

    Analizer.extract_data(event)
    Analizer.clustering_in_towers(merging='off')

    Analizer.Fill1PadEnergy()

    Analizer.clustering_in_towers()

    Analizer.FillNclusters()
    for cluster in range(1, Analizer.get_n_clusters()):
        Analizer.FillClusterDistance(0, cluster)
        Analizer.FillClusterRatio(0, cluster)
        Analizer.FillClusterInverseRatio(0, cluster)
        Analizer.FillClusterDistVsRatio(0, cluster)

    for cluster in range(0, Analizer.get_n_clusters()):
        Analizer.FillClusterEnergy(cluster)
        Analizer.FillClusterPadPos(cluster)
        Analizer.FillClusterNPads(cluster)


output_file = TFile('output_new.root', 'update')

for key in Analizer.h_dict.keys():
    Analizer.h_dict[key].Write()

    # if Analizer.get_energy_ratio(0, 1) > 0.8 and Analizer.get_pad_distance(0, 1) > 10:
    #    Analizer.PlotCheck(event)
    #    print('Plot:', event.apv_evt)

    # if Analizer.get_energy_ratio(0, 2) > 0.8:
    #    Analizer.PlotCheck(event)
    #    print('Plot:', event.apv_evt)

    # if Analizer.get_energy_ratio(0, 3) > 0.7:
    #    Analizer.PlotCheck(event)
    #    print('Plot:', event.apv_evt)


input('Wait')


# ###BE CAREFULL!!! IT MAY MESS UP EXISTING output.root

# tracker_analysis(0)
# tracker_analysis(1)

# energy_fit(0)
# energy_fit(1)
