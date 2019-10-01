from ROOT import TFile, gROOT, TGraphErrors, TH1F, TGraph, nullptr, TF1
import numpy as np

gROOT.SetBatch(0)
gROOT.SetStyle('ATLAS')

# file_data = TFile.Open("../extracted/extracted_data_5gev.root", 'read')
file_data = TFile.Open("../extracted/extracted_data_5gev_merged.root", 'read')
tree_data = file_data.data

output_file = TFile("output.root", "RECREATE")
output_file.cd()


def plot():
    h1 = TH1F('h1', 'Clst1', 8, 0, 8)
    h2 = TH1F('h2', 'Clst2', 8, 0, 8)
    # for event in tree_data:
    tree_data.Draw("tr1_n_hits>>h1", "", "histo")
    h1.SetLineColor(1)

    tree_data.Draw("tr2_n_hits>>h2", "", "histosame")
    h2.SetLineColor(2)

    input("wait")


def electron_identification():
    """I ll write description later"""
    cut_off = np.arange(0.1, 20.1, 0.1)

    n_gen = 0
    n_gen_reco = np.zeros_like(cut_off)
    n_reco = np.zeros_like(cut_off)

    x = np.array([1., 6.])
    y_err = np.array([0.52, 0.52])

    n_events = tree_data.GetEntries()
    for i, event in enumerate(tree_data):
        if i % 1000 == 0:
            if i == 50000:
                break
            print(i, "event out of", n_events)
        # Good event with clear electron generated and detected by calorimeter and trackers
        if not (event.cal_n_clusters > 0
                and event.cal_cluster_energy[0] > 150
                and 130 < event.cal_cluster_y[0] < 150
                and event.tr1_n_hits > 0
                and event.tr2_n_hits > 0):
            continue
        n_gen += 1

        n_electrons = np.zeros_like(cut_off)
        for y3 in event.cal_cluster_y:
            cluster_matched_tracks = np.zeros_like(cut_off)

            for y1 in event.tr1_hit_y:
                for y2 in event.tr2_hit_y:
                    y = np.array([y1, y2])
                    gr_track = TGraphErrors(2, x, y, nullptr, y_err)
                    gr_track.Fit("pol0", "Q")
                    residual = y3 - gr_track.GetFunction("pol0").Eval(24.)

                    for idx, cut in enumerate(cut_off):
                        if -cut < residual < cut:
                            cluster_matched_tracks[idx] += 1

            for idx, cut in enumerate(cut_off):
                if cluster_matched_tracks[idx] > 0:
                    n_electrons[idx] += 1

        for idx, n_e in enumerate(n_electrons):
            if n_e > 0:
                n_reco[idx] += n_e
                n_gen_reco[idx] += 1

    # Get the graph
    eff = n_gen_reco / n_gen * 100
    purity = n_gen_reco / n_reco * 100

    gr = TGraph(len(cut_off), eff, purity)
    gr.GetXaxis().SetTitle("Efficiency, %")
    gr.GetYaxis().SetTitle("Purity, %")
    gr.SetTitle("Tr1 Tr2 fit electron")

    gr.Write("electron_puref")

    return gr


def photon_identification():
    """I ll write description later"""
    cut_off = np.arange(0.1, 20.1, 0.1)

    n_gen = 0
    n_gen_reco = np.zeros_like(cut_off)
    n_reco = np.zeros_like(cut_off)

    x = np.array([1., 6.])
    y_err = np.array([0.52, 0.52])

    n_events = tree_data.GetEntries()
    for i, event in enumerate(tree_data):
        if i % 1000 == 0:
            if i == 50000:
                break
            print(i, "event out of", n_events)
        # Good event with clear electron generated and detected by calorimeter and trackers
        if not (event.cal_n_clusters > 1
                and event.cal_cluster_energy[0] > 150
                and 130 < event.cal_cluster_y[0] < 150
                and 158 < event.cal_cluster_y[1] < 178
                and event.tr1_n_hits > 0
                and event.tr2_n_hits > 0):
            continue
        n_gen += 1

        n_photons = np.zeros_like(cut_off)
        for y3 in event.cal_cluster_y:
            cluster_matched_tracks = np.zeros_like(cut_off)

            for y1 in event.tr1_hit_y:
                for y2 in event.tr2_hit_y:
                    y = np.array([y1, y2])
                    gr_track = TGraphErrors(2, x, y, nullptr, y_err)
                    gr_track.Fit("pol0", "Q")
                    residual = y3 - gr_track.GetFunction("pol0").Eval(24.)

                    for idx, cut in enumerate(cut_off):
                        if -cut < residual < cut:
                            cluster_matched_tracks[idx] += 1

            for idx, cut in enumerate(cut_off):
                if cluster_matched_tracks[idx] == 0:
                    n_photons[idx] += 1

        for idx, n_ph in enumerate(n_photons):
            if n_ph > 0:
                n_reco[idx] += n_ph
                n_gen_reco[idx] += 1

    # Get the graph
    eff = n_gen_reco / n_gen * 100
    purity = n_gen_reco / n_reco * 100

    gr = TGraph(len(cut_off), eff, purity)
    gr.GetXaxis().SetTitle("Efficiency, %")
    gr.GetYaxis().SetTitle("Purity, %")
    gr.SetTitle("Tr1 Tr2 fit photon")
    gr.SetMarkerColor(2)
    gr.SetMarkerStyle(21)

    gr.Write("photon_puref")

    return gr


gr1 = electron_identification()
gr1.Draw("AP")

gr2 = photon_identification()
gr2.Draw("Psame")

input("wait")
