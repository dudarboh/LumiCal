from ROOT import TFile, gROOT, TGraphErrors, TH1F, TGraph, TCanvas, TPad, TF1, gStyle, TF1, nullptr
import numpy as np

gROOT.SetBatch(0)
gROOT.SetStyle('ATLAS')

file_mc_5gev = TFile.Open("../extracted_root_files/extracted_mc_5gev.root", 'read')
tree_mc_5gev = file_mc_5gev.mc

file_mc_4gev = TFile.Open("../extracted_root_files/extracted_mc_4gev.root", 'read')
tree_mc_4gev = file_mc_4gev.mc

file_mc_3gev = TFile.Open("../extracted_root_files/extracted_mc_3gev.root", 'read')
tree_mc_3gev = file_mc_3gev.mc

file_mc_2gev = TFile.Open("../extracted_root_files/extracted_mc_2gev.root", 'read')
tree_mc_2gev = file_mc_2gev.mc

file_mc_1gev = TFile.Open("../extracted_root_files/extracted_mc_1gev.root", 'read')
tree_mc_1gev = file_mc_1gev.mc

# file_5gev = TFile.Open('../extracted_root_files/with_1.4_energy_cut/extracted_5_gev_energy_scan_1.root', 'read')
# tree_5gev = file_5gev.data
file_5gev = TFile.Open("../extracted_root_files/extracted_data_5gev.root", 'read')
tree_5gev = file_5gev.data

file_5gev_no_cd = TFile.Open("../extracted_root_files/extracted_data_nocd_5gev.root", 'read')
tree_5gev_no_cd = file_5gev_no_cd.data

file_4gev = TFile.Open('../extracted_root_files/with_1.4_energy_cut/extracted_4_gev_energy_scan_1.root', 'read')
tree_4gev = file_4gev.data

file_3gev = TFile.Open('../extracted_root_files/with_1.4_energy_cut/extracted_3_gev_energy_scan_1.root', 'read')
tree_3gev = file_3gev.data

file_2gev = TFile.Open('../extracted_root_files/with_1.4_energy_cut/extracted_2_gev_energy_scan_1.root', 'read')
tree_2gev = file_2gev.data

file_1gev = TFile.Open('../extracted_root_files/with_1.4_energy_cut/extracted_1_gev_energy_scan_1.root', 'read')
tree_1gev = file_1gev.data

trees = [tree_1gev, tree_2gev, tree_3gev, tree_4gev, tree_5gev]
# trees = [tree_5gev]

output_file = TFile.Open("./results.root", "RECREATE")
output_file.cd()


def check_alignment(tree):
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
    canvas = TCanvas("check_alignment", "title", 1024, 768)
    canvas.cd()

    file_mc = TFile.Open("../extracted_root_files/extracted_mc_5gev.root", 'read')
    tree_mc = file_mc.mc
    tree_mc.Draw('cal_hit_energy>>h_mc(500, 0, 60)', "cal_hit_layer == 4", "histo")
    h_mc = gROOT.FindObject('h_mc')
    h_mc.Scale(1. / tree_mc.GetEntries())
    h_mc.SetTitle("MC")
    h_mc.SetLineColor(2)
    h_mc.SetLineWidth(3)
    mc_mean = h_mc.GetMean()

    tree_5gev.Draw('cal_hit_energy>>h_data(500, 0, 60)', "cal_hit_layer == 4", "histosame")
    h_data = gROOT.FindObject('h_data')
    h_data.SetLineColor(1)
    h_data.SetTitle("CD")
    h_data.SetLineWidth(3)
    h_data.Scale(1. / tree_5gev.GetEntries())
    data_mean = h_data.GetMean()

    tree_5gev_no_cd.Draw('cal_hit_energy>>h_data_no_cd(500, 0, 60)', "cal_hit_layer == 4", "histosame")
    h_data_no_cd = gROOT.FindObject('h_data_no_cd')
    h_data_no_cd.SetLineColor(4)
    h_data_no_cd.SetTitle("NOCD")
    h_data_no_cd.SetLineWidth(3)
    h_data_no_cd.Scale(1. / tree_5gev_no_cd.GetEntries())
    data_nocd_mean = h_data_no_cd.GetMean()

    canvas.BuildLegend()
    print(mc_mean, "MC mean energy")
    print(data_mean, "data mean energy")
    print(data_nocd_mean, "data nocd mean energy")

    input("stop")


def control_plots():
    canvas = TCanvas("cal_energy_1gev", "title", 1024, 768)
    tree_mc_1gev.Draw('cal_hit_energy>>h_mc_cal_energy(500, 0, 100)', "", "histo")
    h_mc_cal_energy = gROOT.FindObject('h_mc_cal_energy')
    h_mc_cal_energy.Scale(1. / tree_mc_1gev.GetEntries())
    h_mc_cal_energy.SetTitle("MC")
    h_mc_cal_energy.GetXaxis().SetTitle("E_{hit}, MIPs")
    h_mc_cal_energy.GetYaxis().SetTitle("#frac{N_{hit}}{N_{events}}, %")
    h_mc_cal_energy.SetLineColor(9)
    h_mc_cal_energy.SetFillColor(9)

    tree_1gev.Draw('cal_hit_energy>>h_data_cal_energy(500, 0, 100)', "", "histosame")
    h_data_cal_energy = gROOT.FindObject('h_data_cal_energy')
    h_data_cal_energy.Scale(1. / tree_1gev.GetEntries())
    h_data_cal_energy.SetTitle("DATA")
    h_data_cal_energy.SetLineColor(1)
    h_data_cal_energy.SetLineWidth(3)
    canvas.BuildLegend()
    canvas.Write("cal_energy")

    canvas = TCanvas("tr1_energy_1gev", "title", 1024, 768)
    tree_mc_1gev.Draw('tr1_hit_energy>>h_mc_tr1_energy(500, 0, 10)', "", "histo")
    h_mc_tr1_energy = gROOT.FindObject('h_mc_tr1_energy')
    h_mc_tr1_energy.Scale(1. / tree_mc_1gev.GetEntries())
    h_mc_tr1_energy.GetXaxis().SetTitle("E_{hit}, MIPs")
    h_mc_tr1_energy.GetYaxis().SetTitle("#frac{N_{hit}}{N_{events}}, %")
    h_mc_tr1_energy.SetLineColor(9)
    h_mc_tr1_energy.SetFillColor(9)

    tree_1gev.Draw('tr1_hit_energy>>h_data_tr1_energy(500, 0, 10)', "", "histosame")
    h_data_tr1_energy = gROOT.FindObject('h_data_tr1_energy')
    h_data_tr1_energy.Scale(1. / tree_1gev.GetEntries())
    h_data_tr1_energy.SetLineColor(1)
    h_data_tr1_energy.SetLineWidth(3)
    canvas.Write("tr1_energy")

    canvas = TCanvas("tr2_energy_1gev", "title", 1024, 768)
    tree_mc_1gev.Draw('tr2_hit_energy>>h_mc_tr2_energy(500, 0, 10)', "", "histo")
    h_mc_tr2_energy = gROOT.FindObject('h_mc_tr2_energy')
    h_mc_tr2_energy.Scale(1. / tree_mc_1gev.GetEntries())
    h_mc_tr2_energy.GetXaxis().SetTitle("E_{hit}, MIPs")
    h_mc_tr2_energy.GetYaxis().SetTitle("#frac{N_{hit}}{N_{events}}, %")
    h_mc_tr2_energy.SetLineColor(9)
    h_mc_tr2_energy.SetFillColor(9)

    tree_1gev.Draw('tr2_hit_energy>>h_data_tr2_energy(500, 0, 10)', "", "histosame")
    h_data_tr2_energy = gROOT.FindObject('h_data_tr2_energy')
    h_data_tr2_energy.Scale(1. / tree_1gev.GetEntries())
    h_data_tr2_energy.SetLineColor(1)
    h_data_tr2_energy.SetLineWidth(3)
    canvas.Write("tr2_energy")

    canvas = TCanvas("tr2_y_dist_1gev", "title", 1024, 768)
    tree_mc_1gev.Draw('tr2_hit_y - cal_cluster_y[0]>>h_mc_tr2_y_dist(500, -60, 45)', "", "histo")
    h_mc_tr2_y_dist = gROOT.FindObject('h_mc_tr2_y_dist')
    h_mc_tr2_y_dist.Scale(1. / tree_mc_1gev.GetEntries())
    h_mc_tr2_y_dist.GetXaxis().SetTitle("y_{tr2, hit} - y_{main clst}, mm")
    h_mc_tr2_y_dist.GetYaxis().SetTitle("#frac{N_{hit}}{N_{events}}, %")
    h_mc_tr2_y_dist.SetLineColor(9)
    h_mc_tr2_y_dist.SetFillColor(9)

    tree_1gev.Draw('tr2_hit_y - cal_cluster_y[0]>>h_data_tr2_y_dist(500, -60, 45)', "", "histosame")
    h_data_tr2_y_dist = gROOT.FindObject('h_data_tr2_y_dist')
    h_data_tr2_y_dist.Scale(1. / tree_1gev.GetEntries())
    h_data_tr2_y_dist.SetLineColor(1)
    h_data_tr2_y_dist.SetLineWidth(3)
    canvas.Write("tr2_y_dist")

    input("stop")


def std_plot():
    canvas = TCanvas("std_plot", "title", 1024, 768)
    x = [2., 3., 4., 5., 6.]
    y_data_1gev = []
    y_mc_1gev = []
    y_data_2gev = []
    y_mc_2gev = []
    y_data_3gev = []
    y_mc_3gev = []
    y_data_4gev = []
    y_mc_4gev = []
    y_data_5gev = []
    y_mc_5gev = []

    for layer in x:
        tree_mc_1gev.Draw('tr1_hit_energy>>h_mc', "cal_hit_energy*(cal_hit_layer == {})".format(layer))
        h_mc = gROOT.FindObject('h_mc')
        y_mc_1gev.append(h_mc.GetStdDev())

        tree_mc_2gev.Draw('cal_hit_y-cal_cluster_y[0]>>h_mc', "cal_hit_energy*(cal_hit_layer == {})".format(layer))
        h_mc = gROOT.FindObject('h_mc')
        y_mc_2gev.append(h_mc.GetStdDev())

        tree_mc_3gev.Draw('cal_hit_y-cal_cluster_y[0]>>h_mc', "cal_hit_energy*(cal_hit_layer == {})".format(layer))
        h_mc = gROOT.FindObject('h_mc')
        y_mc_3gev.append(h_mc.GetStdDev())

        tree_mc_4gev.Draw('cal_hit_y-cal_cluster_y[0]>>h_mc', "cal_hit_energy*(cal_hit_layer == {})".format(layer))
        h_mc = gROOT.FindObject('h_mc')
        y_mc_4gev.append(h_mc.GetStdDev())

        tree_mc_5gev.Draw('cal_hit_y-cal_cluster_y[0]>>h_mc', "cal_hit_energy*(cal_hit_layer == {})".format(layer))
        h_mc = gROOT.FindObject('h_mc')
        y_mc_5gev.append(h_mc.GetStdDev())

        # Data
        tree_1gev.Draw('cal_hit_y-cal_cluster_y[0]>>h_data', "cal_hit_energy*(cal_hit_layer == {})".format(layer))
        h_data = gROOT.FindObject('h_data')
        y_data_1gev.append(h_data.GetStdDev())

        tree_2gev.Draw('cal_hit_y-cal_cluster_y[0]>>h_data', "cal_hit_energy*(cal_hit_layer == {})".format(layer))
        h_data = gROOT.FindObject('h_data')
        y_data_2gev.append(h_data.GetStdDev())

        tree_3gev.Draw('cal_hit_y-cal_cluster_y[0]>>h_data', "cal_hit_energy*(cal_hit_layer == {})".format(layer))
        h_data = gROOT.FindObject('h_data')
        y_data_3gev.append(h_data.GetStdDev())

        tree_4gev.Draw('cal_hit_y-cal_cluster_y[0]>>h_data', "cal_hit_energy*(cal_hit_layer == {})".format(layer))
        h_data = gROOT.FindObject('h_data')
        y_data_4gev.append(h_data.GetStdDev())

        tree_5gev.Draw('cal_hit_y-cal_cluster_y[0]>>h_data', "cal_hit_energy*(cal_hit_layer == {})".format(layer))
        h_data = gROOT.FindObject('h_data')
        y_data_5gev.append(h_data.GetStdDev())

    gr_data_1gev = TGraphErrors(5, np.array(x), np.array(y_data_1gev), nullptr, nullptr)
    gr_data_1gev.SetTitle("Data 1 GeV")
    gr_data_1gev.SetMarkerStyle(20)
    gr_data_1gev.SetMarkerColor(1)
    gr_data_1gev.Draw("AP")

    gr_data_2gev = TGraphErrors(5, np.array(x), np.array(y_data_2gev), nullptr, nullptr)
    gr_data_2gev.SetTitle("Data 2 GeV")
    gr_data_2gev.SetMarkerStyle(21)
    gr_data_2gev.SetMarkerColor(2)
    gr_data_2gev.Draw("Psame")

    gr_data_3gev = TGraphErrors(5, np.array(x), np.array(y_data_3gev), nullptr, nullptr)
    gr_data_3gev.SetTitle("Data 3 GeV")
    gr_data_3gev.SetMarkerStyle(22)
    gr_data_3gev.SetMarkerColor(8)
    gr_data_3gev.Draw("Psame")

    gr_data_4gev = TGraphErrors(5, np.array(x), np.array(y_data_4gev), nullptr, nullptr)
    gr_data_4gev.SetTitle("Data 4 GeV")
    gr_data_4gev.SetMarkerStyle(23)
    gr_data_4gev.SetMarkerColor(4)
    gr_data_4gev.Draw("Psame")

    gr_data_5gev = TGraphErrors(5, np.array(x), np.array(y_data_5gev), nullptr, nullptr)
    gr_data_5gev.SetTitle("Data 5 GeV")
    gr_data_5gev.SetMarkerStyle(29)
    gr_data_5gev.SetMarkerColor(6)
    gr_data_5gev.Draw("Psame")

    gr_mc_1gev = TGraphErrors(5, np.array(x), np.array(y_mc_1gev), nullptr, nullptr)
    gr_mc_1gev.SetTitle("MC 1 GeV")
    gr_mc_1gev.SetMarkerStyle(1)
    gr_mc_1gev.SetLineColor(1)
    gr_mc_1gev.Draw("Lsame")

    gr_mc_2gev = TGraphErrors(5, np.array(x), np.array(y_mc_2gev), nullptr, nullptr)
    gr_mc_2gev.SetTitle("MC 2 GeV")
    gr_mc_2gev.SetMarkerStyle(1)
    gr_mc_2gev.SetLineColor(2)
    gr_mc_2gev.Draw("Lsame")

    gr_mc_3gev = TGraphErrors(5, np.array(x), np.array(y_mc_3gev), nullptr, nullptr)
    gr_mc_3gev.SetTitle("MC 3 GeV")
    gr_mc_3gev.SetMarkerStyle(1)
    gr_mc_3gev.SetLineColor(8)
    gr_mc_3gev.Draw("Lsame")

    gr_mc_4gev = TGraphErrors(5, np.array(x), np.array(y_mc_4gev), nullptr, nullptr)
    gr_mc_4gev.SetTitle("MC 4 GeV")
    gr_mc_4gev.SetMarkerStyle(1)
    gr_mc_4gev.SetLineColor(4)
    gr_mc_4gev.Draw("Lsame")

    gr_mc_5gev = TGraphErrors(5, np.array(x), np.array(y_mc_5gev), nullptr, nullptr)
    gr_mc_5gev.SetTitle("MC 5 GeV")
    gr_mc_5gev.SetMarkerStyle(1)
    gr_mc_5gev.SetLineColor(6)
    gr_mc_5gev.Draw("Lsame")

    canvas.BuildLegend()
    canvas.Write("std_plot")
    input("stop")


def n_particles_in_cell():
    canvas = TCanvas("ratio_particles_in_cell_1gev", "title", 1024, 768)
    # canvas.Divide(3, 1)
    # canvas.cd(1)
    tree_mc_1gev.Draw('cal_hit_pad:cal_hit_layer>>h_tot(5, 2., 7., 44, 20., 64.)', "cal_hit_n_bs_particles/(cal_hit_n_dir_particles+cal_hit_n_bs_particles)", "colztext")
    h_tot = gROOT.FindObject("h_tot")
    h_tot.Scale(1. / tree_mc_1gev.GetEntries())
    h_tot.SetTitle("Ratio of BS hits in cells")
    h_tot.GetXaxis().SetTitle("layer number")
    h_tot.GetYaxis().SetTitle("pad number")

    h_proj_1gev = h_tot.ProjectionY("h_proj_1gev", 0, 1)
    # canvas.cd(2)
    # tree_mc_1gev.Draw('cal_hit_pad:cal_hit_layer>>h_dir(5, 2., 7., 44, 20., 64.)', "cal_hit_n_dir_particles", "colztext")
    # h_dir = gROOT.FindObject("h_dir")
    # h_dir.Scale(1. / tree_mc_1gev.GetEntries())
    # h_dir.SetTitle("DIRECT hits in cells")
    # h_dir.GetXaxis().SetTitle("layer number")
    # h_dir.GetYaxis().SetTitle("pad number")
    # canvas.cd(3)
    # tree_mc_1gev.Draw('cal_hit_pad:cal_hit_layer>>h_bs(5, 2., 7., 44, 20., 64.)', "cal_hit_n_bs_particles", "colztext")
    # h_bs = gROOT.FindObject("h_bs")
    # h_bs.Scale(1. / tree_mc_1gev.GetEntries())
    # h_bs.SetTitle("BS hits in cells")
    # h_bs.GetXaxis().SetTitle("layer number")
    # h_bs.GetYaxis().SetTitle("pad number")
    canvas.Write("ratio_particles_in_cell_1gev")

    canvas = TCanvas("ratio_particles_in_cell_2gev", "title", 1024, 768)
    # canvas.Divide(3, 1)
    # canvas.cd(1)
    tree_mc_2gev.Draw('cal_hit_pad:cal_hit_layer>>h_tot(5, 2., 7., 44, 20., 64.)', "cal_hit_n_bs_particles/(cal_hit_n_dir_particles+cal_hit_n_bs_particles)", "colztext")
    h_tot = gROOT.FindObject("h_tot")
    h_tot.Scale(1. / tree_mc_2gev.GetEntries())
    h_tot.SetTitle("Ratio of BS hits in cells")
    h_tot.GetXaxis().SetTitle("layer number")
    h_tot.GetYaxis().SetTitle("pad number")

    h_proj_2gev = h_tot.ProjectionY("h_proj_2gev", 0, 1)

    # canvas.cd(2)
    # tree_mc_2gev.Draw('cal_hit_pad:cal_hit_layer>>h_dir(5, 2., 7., 44, 20., 64.)', "cal_hit_n_dir_particles", "colztext")
    # h_dir = gROOT.FindObject("h_dir")
    # h_dir.Scale(1. / tree_mc_2gev.GetEntries())
    # h_dir.SetTitle("DIRECT hits in cells")
    # h_dir.GetXaxis().SetTitle("layer number")
    # h_dir.GetYaxis().SetTitle("pad number")
    # canvas.cd(3)
    # tree_mc_2gev.Draw('cal_hit_pad:cal_hit_layer>>h_bs(5, 2., 7., 44, 20., 64.)', "cal_hit_n_bs_particles", "colztext")
    # h_bs = gROOT.FindObject("h_bs")
    # h_bs.Scale(1. / tree_mc_2gev.GetEntries())
    # h_bs.SetTitle("BS hits in cells")
    # h_bs.GetXaxis().SetTitle("layer number")
    # h_bs.GetYaxis().SetTitle("pad number")
    canvas.Write("ratio_particles_in_cell_2gev")

    canvas = TCanvas("ratio_particles_in_cell_3gev", "title", 1024, 768)
    # canvas.Divide(3, 1)
    # canvas.cd(1)
    tree_mc_3gev.Draw('cal_hit_pad:cal_hit_layer>>h_tot(5, 2., 7., 44, 20., 64.)', "cal_hit_n_bs_particles/(cal_hit_n_dir_particles+cal_hit_n_bs_particles)", "colztext")
    h_tot = gROOT.FindObject("h_tot")
    h_tot.Scale(1. / tree_mc_3gev.GetEntries())
    h_tot.SetTitle("Ratio of BS hits in cells")
    h_tot.GetXaxis().SetTitle("layer number")
    h_tot.GetYaxis().SetTitle("pad number")

    h_proj_3gev = h_tot.ProjectionY("h_proj_3gev", 0, 1)
    # canvas.cd(2)
    # tree_mc_3gev.Draw('cal_hit_pad:cal_hit_layer>>h_dir(5, 2., 7., 44, 20., 64.)', "cal_hit_n_dir_particles", "colztext")
    # h_dir = gROOT.FindObject("h_dir")
    # h_dir.Scale(1. / tree_mc_3gev.GetEntries())
    # h_dir.SetTitle("DIRECT hits in cells")
    # h_dir.GetXaxis().SetTitle("layer number")
    # h_dir.GetYaxis().SetTitle("pad number")
    # canvas.cd(3)
    # tree_mc_3gev.Draw('cal_hit_pad:cal_hit_layer>>h_bs(5, 2., 7., 44, 20., 64.)', "cal_hit_n_bs_particles", "colztext")
    # h_bs = gROOT.FindObject("h_bs")
    # h_bs.Scale(1. / tree_mc_3gev.GetEntries())
    # h_bs.SetTitle("BS hits in cells")
    # h_bs.GetXaxis().SetTitle("layer number")
    # h_bs.GetYaxis().SetTitle("pad number")
    canvas.Write("ratio_particles_in_cell_3gev")

    canvas = TCanvas("ratio_particles_in_cell_4gev", "title", 1024, 768)
    # canvas.Divide(3, 1)
    # canvas.cd(1)
    tree_mc_4gev.Draw('cal_hit_pad:cal_hit_layer>>h_tot(5, 2., 7., 44, 20., 64.)', "cal_hit_n_bs_particles/(cal_hit_n_dir_particles+cal_hit_n_bs_particles)", "colztext")
    h_tot = gROOT.FindObject("h_tot")
    h_tot.Scale(1. / tree_mc_4gev.GetEntries())
    h_tot.SetTitle("Ratio of BS hits in cells")
    h_tot.GetXaxis().SetTitle("layer number")
    h_tot.GetYaxis().SetTitle("pad number")

    h_proj_4gev = h_tot.ProjectionY("h_proj_4gev", 0, 1)
    # canvas.cd(2)
    # tree_mc_4gev.Draw('cal_hit_pad:cal_hit_layer>>h_dir(5, 2., 7., 44, 20., 64.)', "cal_hit_n_dir_particles", "colztext")
    # h_dir = gROOT.FindObject("h_dir")
    # h_dir.Scale(1. / tree_mc_4gev.GetEntries())
    # h_dir.SetTitle("DIRECT hits in cells")
    # h_dir.GetXaxis().SetTitle("layer number")
    # h_dir.GetYaxis().SetTitle("pad number")
    # canvas.cd(3)
    # tree_mc_4gev.Draw('cal_hit_pad:cal_hit_layer>>h_bs(5, 2., 7., 44, 20., 64.)', "cal_hit_n_bs_particles", "colztext")
    # h_bs = gROOT.FindObject("h_bs")
    # h_bs.Scale(1. / tree_mc_4gev.GetEntries())
    # h_bs.SetTitle("BS hits in cells")
    # h_bs.GetXaxis().SetTitle("layer number")
    # h_bs.GetYaxis().SetTitle("pad number")
    canvas.Write("ratio_particles_in_cell_4gev")

    canvas = TCanvas("ratio_particles_in_cell_5gev", "title", 1024, 768)
    # canvas.Divide(3, 1)
    # canvas.cd(1)
    tree_mc_5gev.Draw('cal_hit_pad:cal_hit_layer>>h_tot(5, 2., 7., 44, 20., 64.)', "cal_hit_n_bs_particles/(cal_hit_n_dir_particles+cal_hit_n_bs_particles)", "colztext")
    h_tot = gROOT.FindObject("h_tot")
    h_tot.Scale(1. / tree_mc_5gev.GetEntries())
    h_tot.SetTitle("Ratio of BS hits in cells")
    h_tot.GetXaxis().SetTitle("layer number")
    h_tot.GetYaxis().SetTitle("pad number")

    h_proj_5gev = h_tot.ProjectionY("h_proj_5gev", 0, 1)
    # canvas.cd(2)
    # tree_mc_5gev.Draw('cal_hit_pad:cal_hit_layer>>h_dir(5, 2., 7., 44, 20., 64.)', "cal_hit_n_dir_particles", "colztext")
    # h_dir = gROOT.FindObject("h_dir")
    # h_dir.Scale(1. / tree_mc_5gev.GetEntries())
    # h_dir.SetTitle("DIRECT hits in cells")
    # h_dir.GetXaxis().SetTitle("layer number")
    # h_dir.GetYaxis().SetTitle("pad number")
    # canvas.cd(3)
    # tree_mc_5gev.Draw('cal_hit_pad:cal_hit_layer>>h_bs(5, 2., 7., 44, 20., 64.)', "cal_hit_n_bs_particles", "colztext")
    # h_bs = gROOT.FindObject("h_bs")
    # h_bs.Scale(1. / tree_mc_5gev.GetEntries())
    # h_bs.SetTitle("BS hits in cells")
    # h_bs.GetXaxis().SetTitle("layer number")
    # h_bs.GetYaxis().SetTitle("pad number")
    canvas.Write("ratio_particles_in_cell_5gev")

    h_proj_5gev.Draw("histo")
    h_proj_4gev.Draw("histosame")
    h_proj_3gev.Draw("histosame")
    h_proj_2gev.Draw("histosame")
    h_proj_1gev.Draw("histosame")
    canvas.Write("Projections")

    input("stop")


def avr_tr_energy():
    canvas = TCanvas("avr_tr2_energy", "title", 1024, 768)
    x = [1., 2., 3., 4., 5.]
    y_data = []
    y_mc = []

    tree_mc_1gev.Draw('Sum$(tr2_hit_energy)>>h_mc')
    h_mc = gROOT.FindObject('h_mc')
    y_mc.append(h_mc.GetMean())

    tree_mc_2gev.Draw('Sum$(tr2_hit_energy)>>h_mc')
    h_mc = gROOT.FindObject('h_mc')
    y_mc.append(h_mc.GetMean())

    tree_mc_3gev.Draw('Sum$(tr2_hit_energy)>>h_mc')
    h_mc = gROOT.FindObject('h_mc')
    y_mc.append(h_mc.GetMean())

    tree_mc_4gev.Draw('Sum$(tr2_hit_energy)>>h_mc')
    h_mc = gROOT.FindObject('h_mc')
    y_mc.append(h_mc.GetMean())

    tree_mc_5gev.Draw('Sum$(tr2_hit_energy)>>h_mc')
    h_mc = gROOT.FindObject('h_mc')
    y_mc.append(h_mc.GetMean())

    # Data
    tree_1gev.Draw('Sum$(tr2_hit_energy)>>h_data')
    h_data = gROOT.FindObject('h_data')
    y_data.append(h_data.GetMean())

    tree_2gev.Draw('Sum$(tr2_hit_energy)>>h_data')
    h_data = gROOT.FindObject('h_data')
    y_data.append(h_data.GetMean())

    tree_3gev.Draw('Sum$(tr2_hit_energy)>>h_data')
    h_data = gROOT.FindObject('h_data')
    y_data.append(h_data.GetMean())

    tree_4gev.Draw('Sum$(tr2_hit_energy)>>h_data')
    h_data = gROOT.FindObject('h_data')
    y_data.append(h_data.GetMean())

    tree_5gev.Draw('Sum$(tr2_hit_energy)>>h_data')
    h_data = gROOT.FindObject('h_data')
    y_data.append(h_data.GetMean())

    gr_data = TGraphErrors(5, np.array(x), np.array(y_data), nullptr, nullptr)
    gr_data.SetTitle("Data")
    gr_data.SetMarkerStyle(20)
    gr_data.SetMarkerColor(1)
    gr_data.Draw("AP")

    gr_mc = TGraphErrors(5, np.array(x), np.array(y_mc), nullptr, nullptr)
    gr_mc.SetTitle("MC")
    gr_mc.SetMarkerStyle(1)
    gr_mc.SetLineColor(2)
    gr_mc.Draw("Lsame")

    canvas.BuildLegend()
    canvas.Write("avr_tr2_energy")
    input("stop")


def bs_plot(tree):
    canvas.cd()

    tree.Draw("tr2_hit_energy>>h_all(100, 0., 2.5)", "", "HISTO")
    h_all = gROOT.FindObject("h_all")
    h_all.SetLineColor(1)
    h_all.SetLineWidth(3)
    h_all.SetTitle("All hits")
    h_all.GetYaxis().SetTitle("Number of Hits")
    h_all.GetXaxis().SetTitle("change")

    # tree.Draw("hit_energy>>h_bs(100, 0., 2.5)", "(hit_layer == 1) && (hit_bs == 1)", "HISTOsame")
    # h_bs = gROOT.FindObject("h_bs")
    # h_bs.SetLineColor(2)
    # h_bs.SetLineWidth(3)
    # h_bs.SetTitle("Back-scattered hits")

    # tree.Draw("hit_energy>>h_dir(100, 0., 2.5)", "(hit_layer == 1) && (hit_bs == 0)", "HISTOsame")
    # h_dir = gROOT.FindObject("h_dir")
    # h_dir.SetLineColor(4)
    # h_dir.SetLineWidth(3)
    # h_dir.SetTitle("Direct hits")

    canvas.BuildLegend()
    input("stop")


def bs_numbers(tree):
    tree.Draw(">>list1", "Sum$(hit_layer == 1) > 1")
    list1 = gROOT.FindObject("list1")
    print(100 * list1.GetN() / tree.GetEntries(), " - events with more than one hit in tr2")

    tree.Draw(">>list2", "Sum$((hit_layer == 1) && (hit_bs == 1)) > 0")
    list2 = gROOT.FindObject("list2")
    print(100 * list2.GetN() / tree.GetEntries(), " - events with at least 1 BS hit in tr2")

    tree.Draw(">>list3", "Sum$((hit_layer == 1) && (hit_bs == 0)) > 1")
    list3 = gROOT.FindObject("list3")
    print(100 * list3.GetN() / tree.GetEntries(), " - events with at least 2 direct hits in tr2")

    input("stop")


def bs_selection_plot(tree):
    canvas.cd()

    tree.Draw("tr2_hit_sector>>h_all(4, 0, 4)", "", "HISTO")
    h_all = gROOT.FindObject("h_all")
    h_all.SetLineColor(1)
    h_all.SetLineWidth(3)
    h_all.SetTitle("All hits")
    h_all.GetYaxis().SetTitle("Number of Hits")
    h_all.GetXaxis().SetTitle("change")

    tree.Draw("tr2_hit_sector>>h_bs(4, 0, 4)", "(tr2_hit_bs == 1)", "samehisto")
    h_bs = gROOT.FindObject("h_bs")
    h_bs.SetLineColor(2)
    h_bs.SetLineWidth(3)
    h_bs.SetTitle("Back-scattered hits")

    tree.Draw("tr2_hit_sector>>h_dir(4, 0, 4)", "(tr2_hit_bs == 0)", "samehisto")
    h_dir = gROOT.FindObject("h_dir")
    h_dir.SetLineColor(4)
    h_dir.SetLineWidth(3)
    h_dir.SetTitle("Direct hits")

    canvas.BuildLegend()
    input("stop")


def bs_selection_numbers(tree):

    tree.Draw(">>list1", "Sum$(tr2_hit_layer == 1) > 1")
    list1 = gROOT.FindObject("list1")
    print(100 * list1.GetN() / tree.GetEntries(), " - events with more than one hit in tr2")

    tree.Draw(">>list2", "Sum$((tr2_hit_layer == 1) && (tr2_hit_bs == 1)) > 0")
    list2 = gROOT.FindObject("list2")
    print(100 * list2.GetN() / tree.GetEntries(), " - events with at least 1 BS hit in tr2")

    tree.Draw(">>list3", "Sum$((tr2_hit_layer == 1) && (tr2_hit_bs == 0)) > 1")
    list3 = gROOT.FindObject("list3")
    print(100 * list3.GetN() / tree.GetEntries(), " - events with at least 2 direct hits in tr2")

    input("stop")


def data_plot():
    canvas.cd()

    tree_mc.Draw("cal_hit_energy>>h_mc(300, 0., 40.)", "cal_hit_layer == 2", "histo")
    h_mc = gROOT.FindObject("h_mc")
    h_mc.Scale(1. / tree_mc.GetEntries())
    h_mc.SetLineColor(9)
    h_mc.SetFillColor(9)
    h_mc.SetTitle('Monte Carlo 5 GeV')
    h_mc.GetYaxis().SetTitle("Number of Hits")
    h_mc.GetXaxis().SetTitle("change")
    # podgonian_for_data_to_mc = 1.0840826105153776

    for i, tree in enumerate(trees):
        tree.Draw('cal_hit_energy>>histo{}(300, 0., 40.)'.format(i + 1), "cal_hit_layer == 2", "histosame")
        histo = gROOT.FindObject('histo{}'.format(i + 1))
        histo.Scale(1. / tree.GetEntries())
        histo.SetLineWidth(3)
        histo.SetLineColor(i + 1)
        histo.SetTitle('Beam energy {} GeV'.format(i + 1))

    # func = TF1("func", "[0]*(1 + TMath::Erf((x-[1])/[2]))", 0.1, 0.65)
    # histo.Fit("func", "R")
    # func.Draw("Lsame")
    # func.SetLineColor(2)

    canvas.BuildLegend()
    input("stop")


def data_numbers():
    canvas.cd()

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


def noise_ana():
    canvas.cd()

    tree_mc.Draw("cal_hit_energy>>h2(500, 0., 10.)", "", "histo")
    h2 = gROOT.FindObject("h2")
    h2.Scale(1. / tree_mc.GetEntries())
    h2.SetLineColor(2)
    h2.SetLineWidth(1)
    h2.SetTitle("MC")
    h2.GetYaxis().SetTitle("Number of Hits")
    h2.GetXaxis().SetTitle("change")

    tree_5gev.Draw("cal_hit_energy>>h1(500, 0., 10.)", "", "histosame")
    h1 = gROOT.FindObject("h1")
    h1.Scale(1. / tree_5gev.GetEntries())
    h1.SetLineColor(1)
    h1.SetLineWidth(1)
    h1.SetTitle("DATA")
    h1.GetYaxis().SetTitle("Number of Hits")
    h1.GetXaxis().SetTitle("change")

    # h1.Divide(h2)
    # h1.Draw("histo")
    # func = TF1("func", "[0]*(1 + TMath::Erf((x-[1])/[2]))", 0., 1.)
    # # func.SetParameter(0, 0.5)
    # h1.Fit("func", "R")
    # func.SetParameter(0, 0.5)
    # func.Draw("same")
    # func.SetLineColor(2)


    canvas.BuildLegend()
    input("stop")

# control_plots()
# std_plot()
# n_particles_in_cell()


# avr_tr_energy()
check_calibration()