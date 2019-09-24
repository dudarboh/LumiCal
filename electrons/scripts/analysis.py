from ROOT import TFile, gROOT, TGraphErrors, TGraph, TF1, gStyle, TF1, nullptr, TH2F
from ROOT import TCanvas, TColor, TGaxis, TH1F, TPad
from ROOT import kBlack, kBlue, kRed

import numpy as np

gROOT.SetStyle('ATLAS')
gROOT.SetBatch(0)

# Data trees
f_data_5gev = TFile.Open("../extracted_root_files/extracted_data_5gev.root", 'read')
f_data_4gev = TFile.Open('../extracted_root_files/extracted_data_4gev.root', 'read')
f_data_3gev = TFile.Open('../extracted_root_files/extracted_data_3gev.root', 'read')
f_data_2gev = TFile.Open('../extracted_root_files/extracted_data_2gev.root', 'read')
f_data_1gev = TFile.Open('../extracted_root_files/extracted_data_1gev.root', 'read')
t_data_5gev = f_data_5gev.data
t_data_4gev = f_data_4gev.data
t_data_3gev = f_data_3gev.data
t_data_2gev = f_data_2gev.data
t_data_1gev = f_data_1gev.data

# MC trees
f_mc_5gev = TFile.Open("../extracted_root_files/extracted_lucas_5gev.root", 'read')
f_mc_4gev = TFile.Open("../extracted_root_files/extracted_lucas_4gev.root", 'read')
f_mc_3gev = TFile.Open("../extracted_root_files/extracted_lucas_3gev.root", 'read')
f_mc_2gev = TFile.Open("../extracted_root_files/extracted_lucas_2gev.root", 'read')
f_mc_1gev = TFile.Open("../extracted_root_files/extracted_lucas_1gev.root", 'read')
t_mc_5gev = f_mc_5gev.mc
t_mc_4gev = f_mc_4gev.mc
t_mc_3gev = f_mc_3gev.mc
t_mc_2gev = f_mc_2gev.mc
t_mc_1gev = f_mc_1gev.mc

output_file = TFile.Open("./results.root", "RECREATE")
output_file.cd()


def func():
    gStyle.SetMarkerStyle(1)

    canvas = TCanvas("name", "title", 1024, 768)
    canvas.cd()
    x = [1., 2., 3., 4., 5.]
    y_data = []
    y_mc = []
    y_err = []

    t_mc_1gev.Draw("Sum$(tr2_hit_n_bs_particles)>>h1", "Sum$(tr2_hit_n_bs_particles > 0 && tr2_hit_n_dir_particles == 0) > 0")
    h1 = gROOT.FindObject("h1")
    y_mc.append(h1.GetMean())
    y_err.append(np.sqrt(h1.GetMean()/h1.GetEntries()))

    t_mc_2gev.Draw("Sum$(tr2_hit_n_bs_particles)>>h2", "Sum$(tr2_hit_n_bs_particles > 0 && tr2_hit_n_dir_particles == 0) > 0")
    h2 = gROOT.FindObject("h2")
    y_mc.append(h2.GetMean())
    y_err.append(np.sqrt(h2.GetMean()/h2.GetEntries()))

    t_mc_3gev.Draw("Sum$(tr2_hit_n_bs_particles)>>h3", "Sum$(tr2_hit_n_bs_particles > 0 && tr2_hit_n_dir_particles == 0) > 0")
    h3 = gROOT.FindObject("h3")
    y_mc.append(h3.GetMean())
    y_err.append(np.sqrt(h3.GetMean()/h3.GetEntries()))

    t_mc_4gev.Draw("Sum$(tr2_hit_n_bs_particles)>>h4", "Sum$(tr2_hit_n_bs_particles > 0 && tr2_hit_n_dir_particles == 0) > 0")
    h4 = gROOT.FindObject("h4")
    y_mc.append(h4.GetMean())
    y_err.append(np.sqrt(h4.GetMean()/h4.GetEntries()))

    t_mc_5gev.Draw("Sum$(tr2_hit_n_bs_particles)>>h5", "Sum$(tr2_hit_n_bs_particles > 0 && tr2_hit_n_dir_particles == 0) > 0")
    h5 = gROOT.FindObject("h5")
    y_mc.append(h5.GetMean())
    y_err.append(np.sqrt(h5.GetMean()/h5.GetEntries()))

    gr_mc = TGraphErrors(5, np.array(x), np.array(y_mc), nullptr, np.array(y_err))
    gr_mc.SetTitle("MC")
    gr_mc.SetMarkerStyle(20)
    gr_mc.SetLineColor(2)
    gr_mc.Draw("APE")


    canvas.BuildLegend()

    input('wait')


def bs_numbers():
    gStyle.SetMarkerStyle(1)

    canvas = TCanvas("name", "title", 1024, 768)
    canvas.cd()
    x = [1., 2., 3., 4., 5.]
    y_mc_bs = []
    y_mc_dir = []

    t_mc_1gev.Draw(">>h1_mc", "Sum$(tr2_hit_n_bs_particles>0) > 0")
    h1_mc = gROOT.FindObject("h1_mc")
    y_mc_bs.append(h1_mc.GetN() / t_mc_1gev.GetEntries() * 100)

    t_mc_2gev.Draw(">>h2_mc", "Sum$(tr2_hit_n_bs_particles>0) > 0")
    h2_mc = gROOT.FindObject("h2_mc")
    y_mc_bs.append(h2_mc.GetN() / t_mc_2gev.GetEntries() * 100)

    t_mc_3gev.Draw(">>h3_mc", "Sum$(tr2_hit_n_bs_particles>0) > 0")
    h3_mc = gROOT.FindObject("h3_mc")
    y_mc_bs.append(h3_mc.GetN() / t_mc_3gev.GetEntries() * 100)

    t_mc_4gev.Draw(">>h4_mc", "Sum$(tr2_hit_n_bs_particles>0) > 0")
    h4_mc = gROOT.FindObject("h4_mc")
    y_mc_bs.append(h4_mc.GetN() / t_mc_4gev.GetEntries() * 100)

    t_mc_5gev.Draw(">>h5_mc", "Sum$(tr2_hit_n_bs_particles>0) > 0")
    h5_mc = gROOT.FindObject("h5_mc")
    y_mc_bs.append(h5_mc.GetN() / t_mc_5gev.GetEntries() * 100)




    t_mc_1gev.Draw(">>h1_mc", "Sum$(tr2_hit_n_dir_particles>0) > 1")
    h1_mc = gROOT.FindObject("h1_mc")
    y_mc_dir.append(h1_mc.GetN() / t_mc_1gev.GetEntries() * 100)

    t_mc_2gev.Draw(">>h2_mc", "Sum$(tr2_hit_n_dir_particles>0) > 1")
    h2_mc = gROOT.FindObject("h2_mc")
    y_mc_dir.append(h2_mc.GetN() / t_mc_2gev.GetEntries() * 100)

    t_mc_3gev.Draw(">>h3_mc", "Sum$(tr2_hit_n_dir_particles>0) > 1")
    h3_mc = gROOT.FindObject("h3_mc")
    y_mc_dir.append(h3_mc.GetN() / t_mc_3gev.GetEntries() * 100)

    t_mc_4gev.Draw(">>h4_mc", "Sum$(tr2_hit_n_dir_particles>0) > 1")
    h4_mc = gROOT.FindObject("h4_mc")
    y_mc_dir.append(h4_mc.GetN() / t_mc_4gev.GetEntries() * 100)

    t_mc_5gev.Draw(">>h5_mc", "Sum$(tr2_hit_n_dir_particles>0) > 1")
    h5_mc = gROOT.FindObject("h5_mc")
    y_mc_dir.append(h5_mc.GetN() / t_mc_5gev.GetEntries() * 100)

    gr_mc_bs = TGraphErrors(5, np.array(x), np.array(y_mc_bs), nullptr, nullptr)
    gr_mc_bs.SetTitle("MC back-scattered")
    gr_mc_bs.SetMarkerStyle(20)
    gr_mc_bs.SetLineColor(1)
    gr_mc_bs.Draw("APL")

    gr_mc = TGraphErrors(5, np.array(x), np.array(y_mc_dir), nullptr, nullptr)
    gr_mc.SetTitle("MC pre-scattered")
    gr_mc.SetMarkerStyle(22)
    gr_mc.SetLineColor(2)
    gr_mc.Draw("PLsame")


    canvas.BuildLegend()

    input('wait')


def create_histo(h_name, tree, variable, selection=""):
    h1 = TH1F(h_name, ("Data; E_{hit, tr1}, MIP, mm; #frac{N_{hits}}{N_{events}}, %"), 100, 0, 10)
    h1.SetMarkerStyle(1)
    if tree.GetName() == 'data':
        h1.SetLineColor(1)
    else:
        h1.SetLineColor(2)
        h1.SetTitle("MC")
    h1.SetLineWidth(2)
    h1.SetStats(0)
    tree.Draw("{}>>{}".format(variable, h_name), selection)
    h1.Scale(1. / tree.GetEntries())
    return h1


def createRatio(h1, h2):
    h3 = h1.Clone("h3")
    h3.SetLineColor(kBlack)
    h3.SetMarkerStyle(21)
    h3.SetTitle("")
    h3.SetMinimum(0.8)
    h3.SetMaximum(1.35)
    # Set up plot for markers and errors
    # h3.Sumw2()
    h3.SetStats(0)
    h3.Divide(h2)

    # Adjust y-axis settings
    y = h3.GetYaxis()
    y.SetTitle("Ratio")
    y.SetNdivisions(306)
    y.SetTitleSize(30)
    y.SetTitleFont(43)
    y.SetTitleOffset(1.55)
    y.SetLabelFont(43)
    y.SetLabelSize(25)

    # Adjust x-axis settings
    x = h3.GetXaxis()
    x.SetTitleSize(30)
    x.SetTitleFont(43)
    x.SetTitleOffset(2.5)
    x.SetLabelFont(43)
    x.SetLabelSize(25)

    return h3


def createCanvasPads():
    c = TCanvas("c", "canvas", 800, 800)
    # Upper histogram plot is pad1
    pad1 = TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0)  # joins upper and lower plot
    pad1.SetGridx()
    pad1.Draw()
    # Lower ratio plot is pad2
    c.cd()  # returns to main canvas before defining pad2
    pad2 = TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    pad2.SetTopMargin(0)  # joins upper and lower plot
    pad2.SetBottomMargin(0.2)
    pad2.SetGridx()
    pad2.SetGridy()
    pad2.Draw()
    return c, pad1, pad2


def ratioplot():
    # create required parts
    h1 = create_histo("h1", t_data_5gev, "tr1_hit_energy", "")
    h2 = create_histo("h2", t_mc_5gev, "tr1_hit_energy", "")
    h3 = createRatio(h1, h2)
    c, pad1, pad2 = createCanvasPads()

    # draw everything
    pad1.cd()
    h1.Draw("histo")
    h2.Draw("histosame")

    pad2.cd()
    h3.Draw("ep")

    # To hold window open when running from command line
    input('i am finished')


# func()
bs_numbers()
