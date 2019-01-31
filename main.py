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

from ROOT import TFile, gROOT
import time

from detector import Calorimeter, Tracker
from signal import extract_signal
from cluster import clustering_in_towers
from histos import h_dict

# To measure time of execution
start_time = time.time()

# Dont draw pictures in the end
gROOT.SetBatch(1)


# Extraction
def coincidence_analysis(self):
    shower_pad, shower_sector = self.clusters_calorimeter[0].pad, self.clusters_calorimeter[0].sector

    for tr1_cluster in self.clusters_tracker1:
        if shower_pad-1 < tr1_cluster.pad < shower_pad+1:
            self.n_events_tr1_cal_match += 1
            for tr2_cluster in self.clusters_tracker2:
                if shower_pad-1 < tr2_cluster.pad < shower_pad+1:
                    self.n_events_all_match += 1
                    break
            else:
                return 0
            break
    else:
        return 0
    if len(self.clusters_tracker1) == 1 and len(self.clusters_tracker2) == 1:
        self.n_events_only_shower_match += 1
        return 0
    return 1


def main():

    filename = './trees/run741_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root'
    file = TFile.Open(filename)
    tree = file.apv_reco

    n_events = tree.GetEntries()
    passed_cuts = 0
    have_signal_in_cal = 0

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
        if not clusters_calorimeter:
            continue
        have_signal_in_cal += 1

        calorimeter = Calorimeter(signals_calorimeter, clusters_calorimeter)
        tracker1 = Tracker(signals_tracker1, clusters_tracker1)
        tracker2 = Tracker(signals_tracker2, clusters_tracker2)

        calorimeter.fill_histos(h_dict)
        tracker1.fill_histos(h_dict, 'Tr1')
        tracker2.fill_histos(h_dict, 'Tr2')

        #check_match = Analizer.coincidence_analysis()
        #if check_match == 0:
        #    continue

        # Analize only events with 1 cluster

    # Update/create output file, where to write all histograms
    output_file = TFile('./Analysis/RENAME.root', 'update')

    # Write all histograms to the output root file
    for name in h_dict:
        h_dict[name].Write()

    # Print Hooraay text
    input('Yaay I am finished :3')


# I believe It is faster if main code is written inside main().
main()
