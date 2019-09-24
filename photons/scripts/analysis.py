from ROOT import TFile, gROOT, TGraphErrors, TH1F, TGraph, TCanvas, TPad, TF1, gStyle, TF1, nullptr, TH2F, TPaletteAxis, gPad, gStyle
import numpy as np

gROOT.SetBatch(0)
gROOT.SetStyle('ATLAS')

# file_data = TFile.Open("../extracted/extracted_data_5gev.root", 'read')
file_data = TFile.Open("../extracted/extracted_data_5gev_merged.root", 'read')
tree_data = file_data.data


def plot():
    h1 = TH1F('h1', 'Clst1', 8, 0, 8)
    h2 = TH1F('h2', 'Clst2', 8, 0, 8)
    # for event in tree_data:
    tree_data.Draw("tr1_n_hits>>h1", "", "histo")
    h1.SetLineColor(1)

    tree_data.Draw("tr2_n_hits>>h2", "", "histosame")
    h2.SetLineColor(2)

    input("wait")


def identification():
    x = np.arange(0, 10, 0.1)
    y1 = []
    y2 = []

    tree_data.Draw(">>n_events1", "cal_n_clusters>0")
    n_events1 = gROOT.FindObject("n_events1")
    n_ev1 = n_events1.GetN()

    tree_data.Draw(">>n_events2", "cal_n_clusters>1")
    n_events2 = gROOT.FindObject("n_events2")
    n_ev2 = n_events2.GetN()


    for acc in x:
        print("Obtaining {} point".format(acc))
        tree_data.Draw(">>list1", "(Sum$(abs(cal_cluster_y[0]-tr1_hit_y)<{0}) > 0) && (Sum$(abs(cal_cluster_y[0]-tr2_hit_y)<{0}) > 0)".format(acc))
        # tree_data.Draw(">>list2", "Sum$(abs(cal_cluster_y[0]-tr2_hit_y)<{}) > 0".format(acc))

        l1 = gROOT.FindObject("list1")
        # l2 = gROOT.FindObject("list2")

        y1.append(l1.GetN() / n_ev1 * 100.)
        # y2.append(l2.GetN() / n_ev1 * 100.)

    gr1 = TGraphErrors(len(x), x, np.array(y1), nullptr, nullptr)
    gr1.Draw("APL")
    gr1.SetTitle("Hits in both trackers")
    gr1.GetXaxis().SetTitle("Accepted track distance, mm")
    gr1.GetYaxis().SetTitle("Accepted events, %")

    # gr2 = TGraphErrors(len(x), x, np.array(y2), nullptr, nullptr)
    # gr2.Draw("PLsame")
    # gr2.SetTitle("Tracker 2 hits")
    # gr2.SetMarkerColor(2)
    # gr2.SetLineColor(2)

    input("wait")

# plot()
identification()
