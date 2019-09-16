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


def plot_energy():
    gStyle.SetMarkerStyle(1)

    canvas = TCanvas("name", "title", 1024, 768)
    canvas.cd()

    h0 = TH1F('h0', 'Data', 500, 0., 10.)

    h1 = TH1F('h1', 'MC', 500, 0., 10.)
    h1.SetLineColor(2)

    t_data_5gev.Draw("Sum$(tr1_hit_energy)>>h0", "", "histo")
    h0.Scale(1. / t_data_5gev.GetEntries())

    t_mc_5gev.Draw("Sum$(tr1_hit_energy)>>h1", "", "histosame")
    h1.Scale(1. / t_mc_5gev.GetEntries())

    canvas.BuildLegend()

    input('wait')


def plot_energy2():
    gStyle.SetMarkerStyle(1)

    canvas = TCanvas("name", "title", 1024, 768)
    canvas.cd()

    h0 = TH1F('h0', 'MC total', 200, 0., 10.)

    h1 = TH1F('h1', 'MC back-scattering', 200, 0., 10.)
    h1.SetLineColor(2)

    t_mc_5gev.Draw("tr2_hit_energy>>h0", "", "histo")
    h0.Scale(1. / t_mc_5gev.GetEntries())

    t_mc_5gev.Draw("tr2_hit_energy>>h1", "(tr2_hit_n_bs_particles > 0 && tr2_hit_n_dir_particles == 0)", "histosame")
    h1.Scale(1. / t_mc_5gev.GetEntries())

    h1.Divide(h0)
    h1.Draw("histo")

    canvas.BuildLegend()

    input('wait')


def plot_pad():
    gStyle.SetMarkerStyle(1)

    canvas = TCanvas("name", "title", 1024, 768)
    canvas.cd()

    h0 = TH1F('h0', 'Data', 64, 0, 64)

    h1 = TH1F('h1', 'MC1', 64, 0, 64)
    h1.SetLineColor(2)

    t_data_5gev.Draw("tr1_hit_pad>>h0", "", "histo")
    h0.Scale(1. / t_data_5gev.GetEntries())

    t_mc_5gev.Draw("tr1_hit_pad>>h1", "", "histosame")
    h1.Scale(1. / t_mc_5gev.GetEntries())

    print("Data mean:", h0.GetMean())
    print("MC mean:", h1.GetMean())


    canvas.BuildLegend()

    input('wait')


def plot_y_dist():
    gStyle.SetMarkerStyle(1)

    canvas = TCanvas("name", "title", 1024, 768)
    canvas.cd()

    h0 = TH1F('h0', 'Data', 64, -60, 45)

    h1 = TH1F('h1', 'MC1', 64, -60, 45)
    h1.SetLineColor(2)

    t_data_5gev.Draw("tr1_hit_y-cal_cluster_y[0]>>h0", "", "histo")
    h0.Scale(1. / t_data_5gev.GetEntries())

    t_mc_5gev.Draw("tr1_hit_y-cal_cluster_y[0]>>h1", "(tr1_hit_n_bs_particles > 0 && tr1_hit_n_dir_particles == 0)", "histosame")
    h1.Scale(1. / t_mc_5gev.GetEntries())

    print("Data mean:", h0.GetMean())
    print("MC mean:", h1.GetMean())


    canvas.BuildLegend()

    input('wait')

def plot_y_dist2():
    gStyle.SetMarkerStyle(1)

    canvas = TCanvas("name", "title", 1024, 768)
    canvas.cd()

    h0 = TH1F('h0', 'MC total', 100, -60, 45)
    h1 = TH1F('h1', 'MC back-scattering', 100, -60, 45)

    t_mc_5gev.Draw("tr2_hit_y-cal_cluster_y[0]>>h0", "")
    t_mc_5gev.Draw("tr2_hit_y-cal_cluster_y[0]>>h1", "(tr2_hit_n_bs_particles > 0 && tr2_hit_n_dir_particles == 0)")


    h1.Divide(h0)
    h1.Draw("histo")


    input('wait')


def plot_y_dist_e_dependence():
    gStyle.SetMarkerStyle(1)

    canvas = TCanvas("name", "title", 1024, 768)
    canvas.cd()

    h0_5gev = TH1F('h0_5gev', 'MC total', 20, 0., 10.)
    h1_5gev = TH1F('h1_5gev', 'MC 5 GeV', 20, 0., 10.)

    h0_4gev = TH1F('h0_4gev', 'MC total', 20, 0., 10.)
    h1_4gev = TH1F('h1_4gev', 'MC 4 GeV', 20, 0., 10.)

    h0_3gev = TH1F('h0_3gev', 'MC total', 20, 0., 10.)
    h1_3gev = TH1F('h1_3gev', 'MC 3 GeV', 20, 0., 10.)

    h0_2gev = TH1F('h0_2gev', 'MC total', 20, 0., 10.)
    h1_2gev = TH1F('h1_2gev', 'MC 2 GeV', 20, 0., 10.)

    h0_1gev = TH1F('h0_1gev', 'MC total', 20, 0., 10.)
    h1_1gev = TH1F('h1_1gev', 'MC 1 GeV', 20, 0., 10.)

    t_mc_5gev.Draw("tr2_hit_energy>>h0_5gev", "")
    t_mc_5gev.Draw("tr2_hit_energy>>h1_5gev", "(tr2_hit_n_bs_particles > 0 && tr2_hit_n_dir_particles == 0)")

    t_mc_4gev.Draw("tr2_hit_energy>>h0_4gev", "")
    t_mc_4gev.Draw("tr2_hit_energy>>h1_4gev", "(tr2_hit_n_bs_particles > 0 && tr2_hit_n_dir_particles == 0)")

    t_mc_3gev.Draw("tr2_hit_energy>>h0_3gev", "")
    t_mc_3gev.Draw("tr2_hit_energy>>h1_3gev", "(tr2_hit_n_bs_particles > 0 && tr2_hit_n_dir_particles == 0)")

    t_mc_2gev.Draw("tr2_hit_energy>>h0_2gev", "")
    t_mc_2gev.Draw("tr2_hit_energy>>h1_2gev", "(tr2_hit_n_bs_particles > 0 && tr2_hit_n_dir_particles == 0)")

    t_mc_1gev.Draw("tr2_hit_energy>>h0_1gev", "")
    t_mc_1gev.Draw("tr2_hit_energy>>h1_1gev", "(tr2_hit_n_bs_particles > 0 && tr2_hit_n_dir_particles == 0)")

    h1_5gev.Divide(h0_5gev)
    h1_5gev.Draw("histo")

    h1_4gev.Divide(h0_4gev)
    h1_4gev.Draw("histosame")
    h1_4gev.SetLineColor(2)

    h1_3gev.Divide(h0_3gev)
    h1_3gev.Draw("histosame")
    h1_3gev.SetLineColor(3)

    h1_2gev.Divide(h0_2gev)
    h1_2gev.Draw("histosame")
    h1_2gev.SetLineColor(4)

    h1_1gev.Divide(h0_1gev)
    h1_1gev.Draw("histosame")
    h1_1gev.SetLineColor(5)
    canvas.BuildLegend()

    input('wait')


def bs_numbers():
    gStyle.SetMarkerStyle(1)

    canvas = TCanvas("name", "title", 1024, 768)
    canvas.cd()
    x = [1., 2., 3., 4., 5.]
    y = []

    t_mc_1gev.Draw(">>h1", "Sum$(tr2_hit_n_bs_particles > 0 && tr2_hit_n_dir_particles == 0) > 0")
    h1 = gROOT.FindObject("h1")
    y.append(h1.GetN()/t_mc_1gev.GetEntries()*100)

    t_mc_2gev.Draw(">>h1", "Sum$(tr2_hit_n_bs_particles > 0 && tr2_hit_n_dir_particles == 0) > 0")
    h1 = gROOT.FindObject("h1")
    y.append(h1.GetN()/t_mc_2gev.GetEntries()*100)

    t_mc_3gev.Draw(">>h1", "Sum$(tr2_hit_n_bs_particles > 0 && tr2_hit_n_dir_particles == 0) > 0")
    h1 = gROOT.FindObject("h1")
    y.append(h1.GetN()/t_mc_3gev.GetEntries()*100)

    t_mc_4gev.Draw(">>h1", "Sum$(tr2_hit_n_bs_particles > 0 && tr2_hit_n_dir_particles == 0) > 0")
    h1 = gROOT.FindObject("h1")
    y.append(h1.GetN()/t_mc_4gev.GetEntries()*100)

    t_mc_5gev.Draw(">>h1", "Sum$(tr2_hit_n_bs_particles > 0 && tr2_hit_n_dir_particles == 0) > 0")
    h1 = gROOT.FindObject("h1")
    y.append(h1.GetN()/t_mc_5gev.GetEntries()*100)

    gr_data = TGraphErrors(5, np.array(x), np.array(y), nullptr, nullptr)
    gr_data.SetTitle("MC")
    gr_data.SetMarkerStyle(20)
    gr_data.SetMarkerColor(1)
    gr_data.Draw("AP")

    canvas.BuildLegend()

    input('wait')


def avr_tr_energy():
    canvas = TCanvas("avr_tr1_energy", "title", 1024, 768)
    x = [1., 2., 3., 4., 5.]
    y_data = []
    y_mc = []
    y_mc_wo = []

    t_mc_1gev.Draw('Sum$(tr1_hit_energy)>>h_mc')
    h_mc = gROOT.FindObject('h_mc')
    y_mc.append(h_mc.GetMean())
    t_mc_1gev.Draw('Sum$(tr1_hit_energy*(tr1_hit_n_bs_particles == 0 || tr1_hit_n_dir_particles > 0))>>h_mc')
    h_mc = gROOT.FindObject('h_mc')
    y_mc_wo.append(h_mc.GetMean())

    t_mc_2gev.Draw('Sum$(tr1_hit_energy)>>h_mc')
    h_mc = gROOT.FindObject('h_mc')
    y_mc.append(h_mc.GetMean())
    t_mc_2gev.Draw('Sum$(tr1_hit_energy*(tr1_hit_n_bs_particles == 0 || tr1_hit_n_dir_particles > 0))>>h_mc')
    h_mc = gROOT.FindObject('h_mc')
    y_mc_wo.append(h_mc.GetMean())

    t_mc_3gev.Draw('Sum$(tr1_hit_energy)>>h_mc')
    h_mc = gROOT.FindObject('h_mc')
    y_mc.append(h_mc.GetMean())
    t_mc_3gev.Draw('Sum$(tr1_hit_energy*(tr1_hit_n_bs_particles == 0 || tr1_hit_n_dir_particles > 0))>>h_mc')
    h_mc = gROOT.FindObject('h_mc')
    y_mc_wo.append(h_mc.GetMean())

    t_mc_4gev.Draw('Sum$(tr1_hit_energy)>>h_mc')
    h_mc = gROOT.FindObject('h_mc')
    y_mc.append(h_mc.GetMean())
    t_mc_4gev.Draw('Sum$(tr1_hit_energy*(tr1_hit_n_bs_particles == 0 || tr1_hit_n_dir_particles > 0))>>h_mc')
    h_mc = gROOT.FindObject('h_mc')
    y_mc_wo.append(h_mc.GetMean())

    t_mc_5gev.Draw('Sum$(tr1_hit_energy)>>h_mc')
    h_mc = gROOT.FindObject('h_mc')
    y_mc.append(h_mc.GetMean())
    t_mc_5gev.Draw('Sum$(tr1_hit_energy*(tr1_hit_n_bs_particles == 0 || tr1_hit_n_dir_particles > 0))>>h_mc')
    h_mc = gROOT.FindObject('h_mc')
    y_mc_wo.append(h_mc.GetMean())
    # Data
    t_data_1gev.Draw('Sum$(tr1_hit_energy)>>h_data')
    h_data = gROOT.FindObject('h_data')
    y_data.append(h_data.GetMean())

    t_data_2gev.Draw('Sum$(tr1_hit_energy)>>h_data')
    h_data = gROOT.FindObject('h_data')
    y_data.append(h_data.GetMean())

    t_data_3gev.Draw('Sum$(tr1_hit_energy)>>h_data')
    h_data = gROOT.FindObject('h_data')
    y_data.append(h_data.GetMean())

    t_data_4gev.Draw('Sum$(tr1_hit_energy)>>h_data')
    h_data = gROOT.FindObject('h_data')
    y_data.append(h_data.GetMean())

    t_data_5gev.Draw('Sum$(tr1_hit_energy)>>h_data')
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

    gr_mc_wo = TGraphErrors(5, np.array(x), np.array(y_mc_wo), nullptr, nullptr)
    gr_mc_wo.SetTitle("MC w/o back-scattering")
    gr_mc_wo.SetMarkerStyle(1)
    gr_mc_wo.SetLineColor(3)
    gr_mc_wo.Draw("Lsame")

    canvas.BuildLegend()
    canvas.Write("avr_tr1_energy")
    input("stop")


def data_2d():
    gStyle.SetMarkerStyle(1)
    gStyle.SetPalette(1)
    canvas = TCanvas("name", "title", 1024, 768)

    h1 = TH2F('h1', 'Data', 64, -60, 45, 100, 0., 10.)

    t_data_5gev.Draw("tr2_hit_energy:(tr2_hit_y-cal_cluster_y[0])>>h1")
    h1.Scale(1. / t_data_5gev.GetEntries())
    h1.Draw("colz")
    h1.GetXaxis().SetTitle("y_{hit,tr2}-y_{shower}, mm")
    h1.GetYaxis().SetTitle("E_{hit,tr2}, MIPs")
    h1.SetMaximum(.1)
    h1.SetMinimum(1e-6)

    input('wait')


def mc_2d():
    gStyle.SetMarkerStyle(1)
    gStyle.SetPalette(1)
    canvas = TCanvas("name", "title", 1024, 768)
    h1 = TH2F('h1', 'Ratio', 20, -60, 45, 20, 0., 10.)
    h2 = TH2F('h2', 'MC total', 20, -60, 45, 20, 0., 10.)
    h3 = TH2F('h3', 'MC3', 20, -60, 45, 20, 0., 10.)

    t_mc_5gev.Draw("tr2_hit_energy:(tr2_hit_y-cal_cluster_y[0])>>h1", "tr2_hit_n_bs_particles > 0 && tr2_hit_n_dir_particles == 0")
    t_mc_5gev.Draw("tr2_hit_energy:(tr2_hit_y-cal_cluster_y[0])>>h2")
    h1.Scale(1. / t_mc_5gev.GetEntries())
    h2.Scale(1. / t_mc_5gev.GetEntries())

    h1.Divide(h2)
    h1.Draw("lego2E")
    h1.GetXaxis().SetTitle("y_{hit,tr2}-y_{shower}, mm")
    h1.GetYaxis().SetTitle("E_{hit,tr2}, MIP")
    h1.SetMaximum(1.)
    h1.SetMinimum(1e-6)

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




# plot_y_dist2()
# plot_energy2()
# ratioplot()
# plot_sector()
# plot_pad()
# mc_2d()
# bs_numbers()
# plot_y_dist_e_dependence()
bs_numbers()
# avr_tr_energy()
