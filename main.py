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

from ROOT import TFile, gROOT, TGraphErrors
import time
import numpy as np

from detector import Calorimeter, Tracker
from signal import extract_signal
from cluster import clustering_in_towers
from histos import h_dict

# To measure time of execution
start_time = time.time()

# Dont draw pictures in the end
gROOT.SetBatch(1)


class Efficiency:
    def __init__(self):
        self.x = np.arange(0, 12, 0.05)
        self.tr1_eff = np.zeros_like(self.x)
        self.tr2_eff = np.zeros_like(self.x)
        self.system_eff = np.zeros_like(self.x)
        self.map_of_events = np.zeros((3, 3))


def track_fit(tr1_hit, tr2_hit, cal_hit):
    x = np.array([1., 26.])
    y = np.array([tr1_hit, cal_hit])
    track = TGraphErrors(2, x, y)
    track.Fit('pol1', "Q")
    fit_func = track.GetFunction('pol1')
    # Return residuals
    #print('p0', fit_func.GetParameter(0))
    #print('p1', fit_func.GetParameter(1))

    return (tr2_hit - fit_func.Eval(6.)), fit_func.GetParameter(0), fit_func.GetParameter(1)


def main():

    file = TFile.Open('./trees/run741_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root')
    tree = file.apv_reco
    n_events = tree.GetEntries()

    efficiency = Efficiency()
    passed_cuts = 0
    have_only_1_cluster_in_cal = 0
    have_only_1_cluster_in_tr1 = 0

    for idx, event in enumerate(tree):
        # For debuging. If you need to loop not over all events
        # if idx == 200:
        #    break

        # Printout to see how many events are proceded/left.
        # And how much time spend per event.
        if idx % (500) == 0:
            time_min = (time.time()-start_time) // 60
            time_sec = (time.time()-start_time) % 60
            print('%(idx)i/%(n_events)i events' % locals(), end=' ')
            print('%(time_min)i min' % locals(), end=' ')
            print('%(time_sec)i sec' % locals())

        # Extract data.
        signals = extract_signal(event)
        # If no hits in calorimeter - skip event
        if not signals:
            continue
        signals_tracker1, signals_tracker2, signals_calorimeter = signals

        passed_cuts += 1

        # Do clustering to signals with merging
        clusters_calorimeter = clustering_in_towers(signals_calorimeter, merge='on')
        clusters_tracker1 = clustering_in_towers(signals_tracker1, merge='off')
        clusters_tracker2 = clustering_in_towers(signals_tracker2, merge='off')

        if len(clusters_calorimeter) != 1:
            continue
        have_only_1_cluster_in_cal += 1

        if len(clusters_tracker1) != 1:
            continue
        have_only_1_cluster_in_tr1 += 1

        if len(clusters_tracker2) == 0:
            continue
        calorimeter = Calorimeter(signals_calorimeter, clusters_calorimeter)
        tracker1 = Tracker(signals_tracker1, clusters_tracker1)
        tracker2 = Tracker(signals_tracker2, clusters_tracker2)

        tracker2.track_par0_fit = track_fit(clusters_tracker1[0].pad, clusters_tracker2[0].pad, clusters_calorimeter[0].pad)[1]
        tracker2.track_par1_fit = track_fit(clusters_tracker1[0].pad, clusters_tracker2[0].pad, clusters_calorimeter[0].pad)[2]

        for cluster in range(min(tracker2.n_clusters(), 2)):
            tracker2.residuals[cluster] = track_fit(clusters_tracker1[0].pad, clusters_tracker2[cluster].pad, clusters_calorimeter[0].pad)[0]
        for hit in signals_tracker2:
            tracker2.residual_hits.append(track_fit(clusters_tracker1[0].pad, hit.pad, clusters_calorimeter[0].pad)[0])

        calorimeter.fill_histos(h_dict)
        tracker1.fill_histos(h_dict, 'Tr1')
        tracker2.fill_histos(h_dict, 'Tr2')

    # Update/create output file, where to write all histograms
    output_file = TFile('./Analysis/RENAME.root', 'update')

    # Write all histograms to the output root file
    for name in h_dict:
        h_dict[name].Write()

    print(have_only_1_cluster_in_cal/n_events, 'Events have ONLY 1 signal in cal:', have_only_1_cluster_in_cal)
    print(have_only_1_cluster_in_tr1/n_events, 'Events have ONLY 1 signal in tr1:', have_only_1_cluster_in_tr1)

    # Print Hooraay text
    input('Yaay I am finished :3')


# I believe It is faster if main code is written inside main().
main()
