from ROOT import TFile, gROOT, TGraphErrors, TH1F, TGraph, TCanvas, TPad, TF1, gStyle
import numpy as np


def check_alignment():
    file = TFile.Open("../extracted_root_files/extracted_data_5gev.root", 'read')
    tree = file.data
    x = np.array([0 * 4.5 + 2.25, 5 * 4.5 + 2.25, 23 * 4.5 + 2.25])

    cuts = "tr1_n_clusters == 1 && tr2_n_clusters == 1 && cal_n_clusters == 1"

    tree.Draw('tr1_cluster_y[0]>>h_tr1(64, 80, 195.2)', cuts)
    h_tr1 = gROOT.FindObject('h_tr1')
    y_tr1 = h_tr1.GetMean()

    tree.Draw('tr2_cluster_y[0]>>h_tr2(64, 80, 195.2)', cuts)
    h_tr2 = gROOT.FindObject('h_tr2')
    y_tr2 = h_tr2.GetMean()

    tree.Draw('cal_cluster_y[0]>>h_cal(64, 80, 195.2)', cuts)
    h_cal = gROOT.FindObject('h_cal')
    y_cal = h_cal.GetMean()

    track = TGraphErrors(3, x, np.array([y_tr1, y_tr2, y_cal]))
    track.Fit('pol0', "Q")
    fit_func = track.GetFunction('pol0')
    print("Tr1 needs alignment for:", y_tr1 - fit_func.Eval(0))
    print("Tr2 needs alignment for:", y_tr2 - fit_func.Eval(0))
    print("Cal needs alignment for:", y_cal - fit_func.Eval(0))


def check_calibration():
    canvas = TCanvas("canvas", "title", 1024, 768)
    canvas.cd()

    file_mc = TFile.Open("../extracted_root_files/extracted_mc_5gev.root", 'read')
    tree_mc = file_mc.mc
    tree_mc.Draw('cal_hit_energy>>h_mc(100, 0, 100)', "cal_hit_layer == 4", "histo")
    h_mc = gROOT.FindObject('h_mc')
    h_mc.SetLineColor(2)
    h_mc.SetLineWidth(3)
    mc_mean = h_mc.GetMean()

    file_5gev = TFile.Open('../extracted_root_files/extracted_data_5gev.root', 'read')
    tree_5gev = file_5gev.data

    tree_5gev.Draw('cal_hit_energy/16.8*19.206>>h_data(100, 0, 100)', "cal_hit_layer == 4", "histosame")
    h_data = gROOT.FindObject('h_data')
    h_data.SetLineColor(1)
    h_data.SetLineWidth(3)
    # h_data.Scale(1. / tree_5gev.GetEntries())
    h_mc.Scale(h_data.GetEntries() / h_mc.GetEntries())
    data_mean = h_data.GetMean()

    canvas.SetLogy()
    print(mc_mean, "MC mean energy")
    print(data_mean, "data mean energy")

    input("stop")


def bs_plot():
    canvas = TCanvas("canvas", "title", 1024, 768)
    canvas.cd()

    file_raw_mc = TFile.Open("../mc_root_files/lucas_tb16_5gev.root", 'read')
    tree_raw_mc = file_raw_mc.LumiCal

    tree_raw_mc.Draw("hit_sector>>h_all(4, 0, 4)", "hit_layer == 1", "HISTO")
    h_all = gROOT.FindObject("h_all")
    h_all.SetLineColor(1)
    h_all.SetLineWidth(3)
    h_all.SetTitle("All hits")
    h_all.GetYaxis().SetTitle("Number of Hits")
    h_all.GetXaxis().SetTitle("change")

    tree_raw_mc.Draw("hit_sector>>h_bs(4, 0, 4)", "(hit_layer == 1) && (hit_bs == 1)", "HISTOsame")
    h_bs = gROOT.FindObject("h_bs")
    h_bs.SetLineColor(2)
    h_bs.SetLineWidth(3)
    h_bs.SetTitle("Back-scattered hits")

    tree_raw_mc.Draw("hit_sector>>h_dir(4, 0, 4)", "(hit_layer == 1) && (hit_bs == 0)", "HISTOsame")
    h_dir = gROOT.FindObject("h_dir")
    h_dir.SetLineColor(4)
    h_dir.SetLineWidth(3)
    h_dir.SetTitle("Direct hits")

    canvas.BuildLegend()
    input("stop")


def bs_numbers():
    file_raw_mc = TFile.Open("../mc_root_files/lucas_tb16_5gev.root", 'read')
    tree_raw_mc = file_raw_mc.LumiCal

    tree_raw_mc.Draw(">>list1", "Sum$(hit_layer == 1) > 1")
    list1 = gROOT.FindObject("list1")
    print(100 * list1.GetN() / tree_raw_mc.GetEntries(), " - events with more than one hit in tr2")

    tree_raw_mc.Draw(">>list2", "Sum$((hit_layer == 1) && (hit_bs == 1)) > 0")
    list2 = gROOT.FindObject("list2")
    print(100 * list2.GetN() / tree_raw_mc.GetEntries(), " - events with at least 1 BS hit in tr2")

    tree_raw_mc.Draw(">>list3", "Sum$((hit_layer == 1) && (hit_bs == 0)) > 1")
    list3 = gROOT.FindObject("list3")
    print(100 * list3.GetN() / tree_raw_mc.GetEntries(), " - events with at least 2 direct hits in tr2")

    tree_raw_mc.Draw(">>list4", "Sum$((hit_layer == 1)) > 1 && Sum$((hit_energy/0.885 <0.7)*(hit_layer == 1)*(hit_pad <35)*(hit_sector == 0 || hit_sector == 3)) > 0")
    list4 = gROOT.FindObject("list4")
    print(100 * list4.GetN() / tree_raw_mc.GetEntries(), " - events passes BS cut")

    input("stop")


def bs_selection_plot():
    canvas = TCanvas("canvas", "title", 1024, 768)
    canvas.cd()

    file_mc = TFile.Open("../extracted_root_files/extracted_mc_5gev.root", 'read')
    tree_mc = file_mc.mc

    tree_mc.Draw("tr2_hit_sector>>h_all(4, 0, 4)", "", "HISTO")
    h_all = gROOT.FindObject("h_all")
    h_all.SetLineColor(1)
    h_all.SetLineWidth(3)
    h_all.SetTitle("All hits")
    h_all.GetYaxis().SetTitle("Number of Hits")
    h_all.GetXaxis().SetTitle("change")

    tree_mc.Draw("tr2_hit_sector>>h_bs(4, 0, 4)", "(tr2_hit_bs == 1)", "samehisto")
    h_bs = gROOT.FindObject("h_bs")
    h_bs.SetLineColor(2)
    h_bs.SetLineWidth(3)
    h_bs.SetTitle("Back-scattered hits")

    tree_mc.Draw("tr2_hit_sector>>h_dir(4, 0, 4)", "(tr2_hit_bs == 0)", "samehisto")
    h_dir = gROOT.FindObject("h_dir")
    h_dir.SetLineColor(4)
    h_dir.SetLineWidth(3)
    h_dir.SetTitle("Direct hits")

    canvas.BuildLegend()
    input("stop")


def bs_selection_numbers():
    file_mc = TFile.Open("../extracted_root_files/extracted_mc_5gev.root", 'read')
    tree_mc = file_mc.mc

    tree_mc.Draw(">>list1", "Sum$(tr2_hit_layer == 1) > 1")
    list1 = gROOT.FindObject("list1")
    print(100 * list1.GetN() / tree_mc.GetEntries(), " - events with more than one hit in tr2")

    tree_mc.Draw(">>list2", "Sum$((tr2_hit_layer == 1) && (tr2_hit_bs == 1)) > 0")
    list2 = gROOT.FindObject("list2")
    print(100 * list2.GetN() / tree_mc.GetEntries(), " - events with at least 1 BS hit in tr2")

    tree_mc.Draw(">>list3", "Sum$((tr2_hit_layer == 1) && (tr2_hit_bs == 0)) > 1")
    list3 = gROOT.FindObject("list3")
    print(100 * list3.GetN() / tree_mc.GetEntries(), " - events with at least 2 direct hits in tr2")

    input("stop")


def data_plot():
    canvas = TCanvas("canvas", "title", 1024, 768)
    canvas.cd()

    file_mc = TFile.Open("../extracted_root_files/extracted_mc_5gev.root", 'read')
    tree_mc = file_mc.mc

    file_5gev = TFile.Open('../extracted_root_files/with_1.4_energy_cut/extracted_5_gev_energy_scan_1.root', 'read')
    tree_5gev = file_5gev.data

    file_4gev = TFile.Open('../extracted_root_files/with_1.4_energy_cut/extracted_4_gev_energy_scan_1.root', 'read')
    tree_4gev = file_4gev.data

    file_3gev = TFile.Open('../extracted_root_files/with_1.4_energy_cut/extracted_3_gev_energy_scan_1.root', 'read')
    tree_3gev = file_3gev.data

    file_2gev = TFile.Open('../extracted_root_files/with_1.4_energy_cut/extracted_2_gev_energy_scan_1.root', 'read')
    tree_2gev = file_2gev.data

    file_1gev = TFile.Open('../extracted_root_files/with_1.4_energy_cut/extracted_1_gev_energy_scan_1.root', 'read')
    tree_1gev = file_1gev.data

    trees = [tree_1gev, tree_2gev, tree_3gev, tree_4gev, tree_5gev]

    tree_mc.Draw("tr2_hit_y - cal_cluster_y[0]>>h_mc(100, -60, 45)", "", "HISTO")
    h_mc = gROOT.FindObject("h_mc")
    h_mc.Scale(1. / tree_mc.GetEntries())
    h_mc.SetLineColor(8)
    h_mc.SetFillColor(8)
    h_mc.SetTitle('Monte Carlo 5 GeV')
    h_mc.GetYaxis().SetTitle("Number of Hits")
    h_mc.GetXaxis().SetTitle("change")
    podgonian_for_data_to_mc = 1.0840826105153776

    for i, tree in enumerate(trees):
        tree.Draw('tr2_hit_y - cal_cluster_y[0]>>histo{}(100, -60, 45)'.format(i + 1), "", "histosame")
        histo = gROOT.FindObject('histo{}'.format(i + 1))
        histo.Scale(1. / tree.GetEntries())
        histo.SetLineWidth(3)
        histo.SetLineColor(i + 1)
        histo.SetTitle('Beam energy {} GeV'.format(i + 1))
    canvas.BuildLegend()
    input("stop")


def data_numbers():
    canvas = TCanvas("canvas", "title", 1024, 768)
    canvas.cd()

    file_mc = TFile.Open("../extracted_root_files/extracted_mc_5gev.root", 'read')
    tree_mc = file_mc.mc

    file_5gev = TFile.Open('../extracted_root_files/with_1.4_energy_cut/extracted_5_gev_energy_scan_1.root', 'read')
    tree_5gev = file_5gev.data

    file_4gev = TFile.Open('../extracted_root_files/with_1.4_energy_cut/extracted_4_gev_energy_scan_1.root', 'read')
    tree_4gev = file_4gev.data

    file_3gev = TFile.Open('../extracted_root_files/with_1.4_energy_cut/extracted_3_gev_energy_scan_1.root', 'read')
    tree_3gev = file_3gev.data

    file_2gev = TFile.Open('../extracted_root_files/with_1.4_energy_cut/extracted_2_gev_energy_scan_1.root', 'read')
    tree_2gev = file_2gev.data

    file_1gev = TFile.Open('../extracted_root_files/with_1.4_energy_cut/extracted_1_gev_energy_scan_1.root', 'read')
    tree_1gev = file_1gev.data

    trees = [tree_1gev, tree_2gev, tree_3gev, tree_4gev, tree_5gev]

    tree_mc.Draw(">>h_mc", "Sum$(tr2_hit_layer == 1) > 1")
    h_mc = gROOT.FindObject("h_mc")
    gr_mc = TGraph(1, np.array([5.]), np.array([100 * h_mc.GetN() / tree_mc.GetEntries()]))
    x = np.array([1., 2., 3., 4., 5.])
    y = np.zeros(5)

    for i, tree in enumerate(trees):
        tree.Draw('>>list{}'.format(i + 1), "Sum$(tr2_hit_layer == 1) > 1")
        list = gROOT.FindObject('list{}'.format(i + 1))
        y[i] = 100 * list.GetN() / tree.GetEntries()

    gr_data = TGraph(5, x, y)
    gr_data.Draw("AP")
    gr_data.SetTitle("Data; Beam Energy, GeV; Events with N_{hit, tr2}>1, %")
    gr_data.SetMarkerStyle(20)
    gr_data.SetMarkerColor(4)
    gr_mc.Draw("Psame")
    gr_mc.SetMarkerStyle(20)
    gr_mc.SetMarkerColor(2)
    gr_mc.SetTitle("MC")

    canvas.BuildLegend()
    input("stop")


gROOT.SetBatch(0)
gROOT.SetStyle('ATLAS')

check_alignment()
# bs_plot()
# bs_selection_plot()
# bs_numbers()
# bs_selection_numbers()

# data_plot()
# data_numbers()

# check_calibration()
