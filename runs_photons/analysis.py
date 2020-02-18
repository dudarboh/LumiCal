from ROOT import TFile, gROOT, TGraphErrors, TGraph, TF1, gStyle, TF1, nullptr, TH2F, TMath
from ROOT import TCanvas, TColor, TGaxis, TH1F, TPad
from ROOT import kBlack, kBlue, kRed

import numpy as np

gROOT.SetStyle('ATLAS')
gROOT.SetBatch(0)
gStyle.SetPalette(1)

f_data_5gev = TFile.Open("./extracted_trees/extracted_data_5gev.root", 'read')
f_mc_5gev = TFile.Open("./extracted_trees/extracted_mc_5gev.root", 'read')
t_data_5gev = f_data_5gev.data
t_mc_5gev = f_mc_5gev.mc


def langaufun(x, par):
    invsq2pi = 0.3989422804014  # (2 pi)^(-1 / 2)
    mpshift = -0.22278298  # Landau maximum location

    np = 100  # number of convolution steps
    sc = 5.  # convolution extends to +-sc Gaussian sigmas

    summ = 0

    mpc = par[1] - mpshift * par[0]

    xlow = x[0] - sc * par[3]
    xupp = x[0] + sc * par[3]

    step = (xupp - xlow) / np

    if par[0] == 0:
        par[0] = 1e-6
    if par[3] == 0:
        par[3] = 1e-6
    for i in range(np // 2):
        xx = xlow + (i + 0.5) * step
        fland = TMath.Landau(xx, mpc, par[0]) / par[0]
        summ += fland * TMath.Gaus(x[0], xx, par[3])

        xx = xupp - (i + 0.5) * step
        fland = TMath.Landau(xx, mpc, par[0]) / par[0]
        summ += fland * TMath.Gaus(x[0], xx, par[3])

    return par[2] * step * summ * invsq2pi / par[3]


def createRatio(h1, h2):
    h3 = h1.Clone("h3")
    h3.SetLineColor(kBlack)
    h3.SetMarkerStyle(21)
    h3.SetTitle("")
    # Set up plot for markers and errors
    # h3.Sumw2()
    h3.SetStats(0)
    h3.Divide(h2)
    h3.SetMinimum(0.9*h3.GetMinimum())
    h3.SetMaximum(1.1*h3.GetMaximum())

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
    x.SetTitleSize(23)
    x.SetTitleFont(43)
    x.SetTitleOffset(2.5)
    x.SetLabelFont(43)
    x.SetLabelSize(25)

    return h3


def createCanvasPads(name):
    c = TCanvas(name, name, 800, 800)
    # Upper histogram plot is pad1
    pad1 = TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0)  # joins upper and lower plot
    pad1.SetGridx()
    pad1.SetGridy()
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



def pad_spectra():
    h_total = TH1F("h_total", "total spectra", 50, 0., 3)

    n_events = t_data_5gev.GetEntries()

    h_arr = []
    fits_arr = []
    fits_means = []
    i = 1
    gr = TGraph()

    t_data_5gev.Draw("tr2_hit_energy>>h_total", "", "histo norm")
    fit_total = TF1("fit_total", langaufun, 0.8, 2, 4)
    fit_total.SetParameters(0.07, 1., 1., 0.1)
    fit_total.SetParNames("Width", "MP", "Area", "GSigma")
    fit_total.SetParLimits(0, 0., 0.1)
    fit_total.SetParLimits(1, 0.8, 1.3)
    fit_total.SetParLimits(2, 0., 1.)
    fit_total.SetParLimits(3, 0., 0.5)

    h_total.Fit("fit_total", "R0")
    total_mean = fit_total.GetParameter(1)

    for pad in range(35, 55):
        for sector in range(1, 3):
            h_arr.append(TH1F("h{}".format(sector * 64 + pad), "pad_spectra", 100, 0., 3))

            fits_arr.append(TF1("fit{}".format(sector * 64 + pad), langaufun, 0.8, 2, 4))
            fits_arr[-1].SetParameters(0.07, 1., 0.03, 0.1)
            fits_arr[-1].SetParNames("Width", "MP", "Area", "GSigma")
            fits_arr[-1].SetParLimits(0, 0., 0.1)
            fits_arr[-1].SetParLimits(1, 0.8, 1.3)
            fits_arr[-1].SetParLimits(2, 0., 1.)
            fits_arr[-1].SetParLimits(3, 0., 0.5)
            h_arr[-1].SetLineColor(i)
            i += 1
            if pad == 0 and sector == 0:
                t_data_5gev.Draw("tr2_hit_energy>>h{}".format(sector*64 + pad), "tr2_hit_pad == {} && tr2_hit_sector == {}".format(pad, sector), "histo norm")
            else:
                t_data_5gev.Draw("tr2_hit_energy>>h{}".format(sector*64 + pad), "tr2_hit_pad == {} && tr2_hit_sector == {}".format(pad, sector), "histosame norm")
                print(pad, "pad")

            h_arr[-1].Fit("fit{}".format(sector * 64 + pad), "R0")
            fits_means.append(fits_arr[-1].GetParameter(1))

            gr.SetPoint(i, pad, fits_means[-1]/total_mean)
            input("wait")

    gr.Draw("AP")
    input("wait")


def y_distance():
    # 0 - mixed
    # 1 - primary
    # 2- electrons
    # 3 - gammas
    # 4 - positrons
    # 5 - hadrons
    h0 = TH1F("h0", "Data", 200, -60., 45.)
    h1 = TH1F("h1", "MC Primary", 200, -60., 45.)
    h1.GetXaxis().SetTitle("y_{tr2,hit} - y_{cal,clst}, mm")
    h1.GetYaxis().SetTitle("N_{hits}/N_{events}")
    h2 = TH1F("h2", "MC PS e^{-}", 200, -60., 45.)
    h3 = TH1F("h3", "MC PS photons", 200, -60., 45.)
    h4 = TH1F("h4", "MC mixed", 200, -60., 45.)
    h5 = TH1F("h5", "MC BS e^{-}", 200, -60., 45.)
    h6 = TH1F("h6", "MC BS photons", 200, -60., 45.)
    h7 = TH1F("h7", "MC BS e^{+}", 200, -60., 45.)
    h8 = TH1F("h8", "MC BS hadrons", 200, -60., 45.)
    histos = [h1, h2, h3, h4, h5, h6, h7, h8]

    bs_hadr = "(tr2_particle_pz < 0 && tr2_hit_type == 5)"
    bs_pos = "(tr2_particle_pz < 0 && tr2_hit_type == 4)"
    bs_gamma = "(tr2_particle_pz < 0 && tr2_hit_type == 3)"
    bs_e = "(tr2_particle_pz < 0 && tr2_hit_type == 2)"

    ps_gamma = "(tr2_particle_pz > 0 && tr2_hit_type == 3)"
    ps_e = "(tr2_particle_pz > 0 && tr2_hit_type == 2)"
    ps_prim = "(tr2_particle_pz > 0 && tr2_hit_type == 1)"
    mix = "(tr2_hit_type == 0)"

    t_data_5gev.Draw("tr2_hit_y - cal_cluster_y[0]>>h0")
    h0.Scale(1. / t_data_5gev.GetEntries())
    t_mc_5gev.Draw("tr2_hit_y - cal_cluster_y[0]>>h8", bs_hadr)
    t_mc_5gev.Draw("tr2_hit_y - cal_cluster_y[0]>>h7", bs_hadr + "||" + bs_pos)
    t_mc_5gev.Draw("tr2_hit_y - cal_cluster_y[0]>>h6", bs_hadr + "||" + bs_pos + "||" + bs_gamma)
    t_mc_5gev.Draw("tr2_hit_y - cal_cluster_y[0]>>h5", bs_hadr + "||" + bs_pos + "||" + bs_gamma + "||" + bs_e)
    t_mc_5gev.Draw("tr2_hit_y - cal_cluster_y[0]>>h4", bs_hadr + "||" + bs_pos + "||" + bs_gamma + "||" + bs_e + "||" + mix)
    t_mc_5gev.Draw("tr2_hit_y - cal_cluster_y[0]>>h3", bs_hadr + "||" + bs_pos + "||" + bs_gamma + "||" + bs_e + "||" + mix + "||" + ps_gamma)
    t_mc_5gev.Draw("tr2_hit_y - cal_cluster_y[0]>>h2", bs_hadr + "||" + bs_pos + "||" + bs_gamma + "||" + bs_e + "||" + mix + "||" + ps_gamma + "||" + ps_e)
    t_mc_5gev.Draw("tr2_hit_y - cal_cluster_y[0]>>h1", bs_hadr + "||" + bs_pos + "||" + bs_gamma + "||" + bs_e + "||" + mix + "||" + ps_gamma + "||" + ps_e + "||" + ps_prim)

    for i, h in enumerate(histos):
        h.SetMarkerStyle(0)
        h.SetLineWidth(1)
        h.Scale(1. / t_mc_5gev.GetEntries())

    # Data divide by MC
    h_ratio = createRatio(h0, h1)
    c, pad1, pad2 = createCanvasPads("tr2_y_distance_ph")

    # draw everything
    pad1.cd()

    h1.Draw("histo")
    h1.SetFillColor(8)
    h2.Draw("histosame")
    h2.SetFillColor(9)
    h3.Draw("histosame")
    h3.SetFillColor(28)
    h4.Draw("histosame")
    h4.SetFillColor(30)
    h5.Draw("histosame")
    h5.SetFillColor(38)
    h6.Draw("histosame")
    h6.SetFillColor(46)
    h7.Draw("histosame")
    h7.SetFillColor(36)

    h0.SetMarkerStyle(20)
    h0.Draw("Psame")

    pad2.cd()
    h_ratio.Draw("ehistop")
    h_ratio.GetXaxis().SetTitle(h1.GetXaxis().GetTitle())

    input("wait")


def dep_energy():
    # 0 - mixed
    # 1 - primary
    # 2- electrons
    # 3 - gammas
    # 4 - positrons
    # 5 - hadrons
    h0 = TH1F("h0", "Data", 200, 0., 10.*0.0885)
    h1 = TH1F("h1", "MC Primary", 200, 0., 10.*0.0885)
    h1.GetXaxis().SetTitle("E_{tr1,hit}, MeV")
    h1.GetYaxis().SetTitle("N_{hits}/N_{events}")
    h2 = TH1F("h2", "MC PS e^{-}", 200, 0., 10.*0.0885)
    h3 = TH1F("h3", "MC PS photons", 200, 0., 10.*0.0885)
    h4 = TH1F("h4", "MC mixed", 200, 0., 10.*0.0885)
    h5 = TH1F("h5", "MC BS e^{-}", 200, 0., 10.*0.0885)
    h6 = TH1F("h6", "MC BS photons", 200, 0., 10.*0.0885)
    h7 = TH1F("h7", "MC BS e^{+}", 200, 0., 10.*0.0885)
    h8 = TH1F("h8", "MC BS hadrons", 200, 0., 10.*0.0885)
    histos = [h1, h2, h3, h4, h5, h6, h7, h8]

    bs_hadr = "(tr2_particle_pz < 0 && tr2_hit_type == 5)"
    bs_pos = "(tr2_particle_pz < 0 && tr2_hit_type == 4)"
    bs_gamma = "(tr2_particle_pz < 0 && tr2_hit_type == 3)"
    bs_e = "(tr2_particle_pz < 0 && tr2_hit_type == 2)"

    ps_gamma = "(tr2_particle_pz > 0 && tr2_hit_type == 3)"
    ps_e = "(tr2_particle_pz > 0 && tr2_hit_type == 2)"
    ps_prim = "(tr2_particle_pz > 0 && tr2_hit_type == 1)"
    mix = "(tr2_hit_type == 0)"

    t_data_5gev.Draw("tr2_cluster_energy*0.0885>>h0")
    h0.Scale(1. / t_data_5gev.GetEntries())
    t_mc_5gev.Draw("tr2_hit_energy*0.0885>>h8", bs_hadr)
    t_mc_5gev.Draw("tr2_hit_energy*0.0885>>h7", bs_hadr + "||" + bs_pos)
    t_mc_5gev.Draw("tr2_hit_energy*0.0885>>h6", bs_hadr + "||" + bs_pos + "||" + bs_gamma)
    t_mc_5gev.Draw("tr2_hit_energy*0.0885>>h5", bs_hadr + "||" + bs_pos + "||" + bs_gamma + "||" + bs_e)
    t_mc_5gev.Draw("tr2_hit_energy*0.0885>>h4", bs_hadr + "||" + bs_pos + "||" + bs_gamma + "||" + bs_e + "||" + mix)
    t_mc_5gev.Draw("tr2_hit_energy*0.0885>>h3", bs_hadr + "||" + bs_pos + "||" + bs_gamma + "||" + bs_e + "||" + mix + "||" + ps_gamma)
    t_mc_5gev.Draw("tr2_hit_energy*0.0885>>h2", bs_hadr + "||" + bs_pos + "||" + bs_gamma + "||" + bs_e + "||" + mix + "||" + ps_gamma + "||" + ps_e)
    t_mc_5gev.Draw("tr2_hit_energy*0.0885>>h1", bs_hadr + "||" + bs_pos + "||" + bs_gamma + "||" + bs_e + "||" + mix + "||" + ps_gamma + "||" + ps_e + "||" + ps_prim)

    for i, h in enumerate(histos):
        h.SetMarkerStyle(0)
        h.SetLineWidth(1)
        h.Scale(1. / t_mc_5gev.GetEntries())

    # Data divide by MC
    h_ratio = createRatio(h0, h1)
    c, pad1, pad2 = createCanvasPads("tr2_dep_energy_ph")

    # draw everything
    pad1.cd()

    h1.Draw("histo")
    h1.SetFillColor(8)
    h2.Draw("histosame")
    h2.SetFillColor(9)
    h3.Draw("histosame")
    h3.SetFillColor(28)
    h4.Draw("histosame")
    h4.SetFillColor(30)
    h5.Draw("histosame")
    h5.SetFillColor(38)
    h6.Draw("histosame")
    h6.SetFillColor(46)
    h7.Draw("histosame")
    h7.SetFillColor(36)

    h0.SetMarkerStyle(20)
    h0.Draw("Psame")

    pad2.cd()
    h_ratio.Draw("ehistop")
    h_ratio.GetXaxis().SetTitle(h1.GetXaxis().GetTitle())

    input("wait")


def particle_energy():
    # 0 - mixed
    # 1 - primary
    # 2- electrons
    # 3 - gammas
    # 4 - positrons
    # 5 - hadrons
    h2 = TH1F("h2", "MC PS e^{-}", 200, 0., 10.)
    h3 = TH1F("h3", "MC PS photons", 200, 0., 10.)
    h4 = TH1F("h4", "MC BS e^{-}", 200, 0., 10.)
    h5 = TH1F("h5", "MC BS photons", 200, 0., 10.)
    h6 = TH1F("h6", "MC BS e^{+}", 200, 0., 10.)
    histos = [h2, h3, h4, h5, h6]

    bs_hadr = "(tr2_particle_pz < 0 && tr2_hit_type == 5)"
    bs_pos = "(tr2_particle_pz < 0 && tr2_hit_type == 4)"
    bs_gamma = "(tr2_particle_pz < 0 && tr2_hit_type == 3)"
    bs_e = "(tr2_particle_pz < 0 && tr2_hit_type == 2)"

    ps_gamma = "(tr2_particle_pz > 0 && tr2_hit_type == 3)"
    ps_e = "(tr2_particle_pz > 0 && tr2_hit_type == 2)"

    t_mc_5gev.Draw("tr2_particle_energy>>h6", bs_hadr + "||" + bs_pos)
    t_mc_5gev.Draw("tr2_particle_energy>>h5", bs_hadr + "||" + bs_pos + "||" + bs_gamma)
    t_mc_5gev.Draw("tr2_particle_energy>>h4", bs_hadr + "||" + bs_pos + "||" + bs_gamma + "||" + bs_e)
    t_mc_5gev.Draw("tr2_particle_energy>>h3", bs_hadr + "||" + bs_pos + "||" + bs_gamma + "||" + bs_e + "||" + ps_gamma)
    t_mc_5gev.Draw("tr2_particle_energy>>h2", bs_hadr + "||" + bs_pos + "||" + bs_gamma + "||" + bs_e + "||" + ps_gamma + "||" + ps_e)

    for i, h in enumerate(histos):
        h.SetMarkerStyle(0)
        h.SetLineWidth(1)
        h.Scale(1. / t_mc_5gev.GetEntries())

    h2.Draw("histo")
    h2.GetXaxis().SetTitle("E_{tr2,particle}, MeV")
    h2.GetYaxis().SetTitle("N_{hits}/N_{events}")
    h2.SetFillColor(9)
    h3.Draw("histosame")
    h3.SetFillColor(28)
    h4.Draw("histosame")
    h4.SetFillColor(30)
    h5.Draw("histosame")
    h5.SetFillColor(38)
    h6.Draw("histosame")
    h6.SetFillColor(46)

    input("wait")


def track_len():
    # 0 - mixed
    # 1 - primary
    # 2- electrons
    # 3 - gammas
    # 4 - positrons
    # 5 - hadrons
    h1 = TH1F("h1", "MC Primary", 200, 0., 1.)
    h1.GetXaxis().SetTitle("l_{tr1, track}, mm")
    h1.GetYaxis().SetTitle("N_{hits}/N_{events}")
    h2 = TH1F("h2", "MC PS e^{-}", 200, 0., 1.)
    h3 = TH1F("h3", "MC PS photons", 200, 0., 1.)
    h4 = TH1F("h4", "MC mixed", 200, 0., 1.)
    h5 = TH1F("h5", "MC BS e^{-}", 200, 0., 1.)
    h6 = TH1F("h6", "MC BS photons", 200, 0., 1.)
    h7 = TH1F("h7", "MC BS e^{+}", 200, 0., 1.)
    h8 = TH1F("h8", "MC BS hadrons", 200, 0., 1.)
    histos = [h1, h2, h3, h4, h5, h6, h7, h8]

    bs_hadr = "(tr2_particle_pz < 0 && tr2_hit_type == 5)"
    bs_pos = "(tr2_particle_pz < 0 && tr2_hit_type == 4)"
    bs_gamma = "(tr2_particle_pz < 0 && tr2_hit_type == 3)"
    bs_e = "(tr2_particle_pz < 0 && tr2_hit_type == 2)"

    ps_gamma = "(tr2_particle_pz > 0 && tr2_hit_type == 3)"
    ps_e = "(tr2_particle_pz > 0 && tr2_hit_type == 2)"
    ps_prim = "(tr2_particle_pz > 0 && tr2_hit_type == 1)"
    mix = "(tr2_hit_type == 0)"

    t_mc_5gev.Draw("tr2_track_len>>h8", bs_hadr)
    t_mc_5gev.Draw("tr2_track_len>>h7", bs_hadr + "||" + bs_pos)
    t_mc_5gev.Draw("tr2_track_len>>h6", bs_hadr + "||" + bs_pos + "||" + bs_gamma)
    t_mc_5gev.Draw("tr2_track_len>>h5", bs_hadr + "||" + bs_pos + "||" + bs_gamma + "||" + bs_e)
    t_mc_5gev.Draw("tr2_track_len>>h4", bs_hadr + "||" + bs_pos + "||" + bs_gamma + "||" + bs_e + "||" + mix)
    t_mc_5gev.Draw("tr2_track_len>>h3", bs_hadr + "||" + bs_pos + "||" + bs_gamma + "||" + bs_e + "||" + mix + "||" + ps_gamma)
    t_mc_5gev.Draw("tr2_track_len>>h2", bs_hadr + "||" + bs_pos + "||" + bs_gamma + "||" + bs_e + "||" + mix + "||" + ps_gamma + "||" + ps_e)
    t_mc_5gev.Draw("tr2_track_len>>h1", bs_hadr + "||" + bs_pos + "||" + bs_gamma + "||" + bs_e + "||" + mix + "||" + ps_gamma + "||" + ps_e + "||" + ps_prim)

    for i, h in enumerate(histos):
        h.SetMarkerStyle(0)
        h.SetLineWidth(1)
        h.Scale(1. / t_mc_5gev.GetEntries())

    h1.Draw("histo")
    h1.SetFillColor(8)
    h2.Draw("histosame")
    h2.SetFillColor(9)
    h3.Draw("histosame")
    h3.SetFillColor(28)
    h4.Draw("histosame")
    h4.SetFillColor(30)
    h5.Draw("histosame")
    h5.SetFillColor(38)
    h6.Draw("histosame")
    h6.SetFillColor(46)
    h7.Draw("histosame")
    h7.SetFillColor(36)

    input("wait")


def calculate_numbers():
    # 0 - mixed
    # 1 - primary
    # 2 - electrons
    # 3 - gammas
    # 4 - positrons
    # 5 - hadrons

    bs_e = "(tr1_particle_pz < 0 && tr1_hit_type == 2)"
    bs_gamma = "(tr1_particle_pz < 0 && tr1_hit_type == 3)"
    bs_pos = "(tr1_particle_pz < 0 && tr1_hit_type == 4)"
    bs_hadr = "(tr1_particle_pz < 0 && tr1_hit_type == 5)"
    any_bs = "(tr1_particle_pz < 0 && tr1_hit_type != 0)"

    ps_e = "(tr1_particle_pz > 0 && tr1_hit_type == 2)"
    ps_gamma = "(tr1_particle_pz > 0 && tr1_hit_type == 3)"
    ps_pos = "(tr1_particle_pz > 0 && tr1_hit_type == 4)"
    ps_hadr = "(tr1_particle_pz > 0 && tr1_hit_type == 5)"
    any_ps = "(tr1_particle_pz > 0 && tr1_hit_type != 1 && tr1_hit_type != 0)"

    ps_prim = "(tr1_particle_pz > 0 && tr1_hit_type == 1)"
    mix = "(tr1_hit_type == 0)"

    t_mc_5gev.Draw("tr1_hit_type>>h_n_tracks")
    h_n_tracks = gROOT.FindObject("h_n_tracks")
    n_tracks = h_n_tracks.GetEntries()
    print("N Tracks: ", n_tracks)

    t_mc_5gev.Draw("tr1_hit_type>>h12", any_ps)
    t_mc_5gev.Draw("tr1_hit_type>>h11", any_bs)
    t_mc_5gev.Draw("tr1_hit_type>>h10", ps_hadr)
    t_mc_5gev.Draw("tr1_hit_type>>h9", ps_pos)
    t_mc_5gev.Draw("tr1_hit_type>>h8", bs_hadr)
    t_mc_5gev.Draw("tr1_hit_type>>h7", bs_pos)
    t_mc_5gev.Draw("tr1_hit_type>>h6", bs_gamma)
    t_mc_5gev.Draw("tr1_hit_type>>h5", bs_e)
    t_mc_5gev.Draw("tr1_hit_type>>h4", mix)
    t_mc_5gev.Draw("tr1_hit_type>>h3", ps_gamma)
    t_mc_5gev.Draw("tr1_hit_type>>h2", ps_e)
    t_mc_5gev.Draw("tr1_hit_type>>h1", ps_prim)

    h1 = gROOT.FindObject("h1")
    print(h1.GetEntries() / n_tracks * 100, " - events with primary particle")
    h2 = gROOT.FindObject("h2")
    print(h2.GetEntries() / n_tracks * 100, " - events with PS e^{-} signal")
    h3 = gROOT.FindObject("h3")
    print(h3.GetEntries() / n_tracks * 100, " - events with PS photon signal")
    h4 = gROOT.FindObject("h4")
    print(h4.GetEntries() / n_tracks * 100, " - events with mixed partcle signal")
    h5 = gROOT.FindObject("h5")
    print(h5.GetEntries() / n_tracks * 100, " - events with BS e^- signal")
    h6 = gROOT.FindObject("h6")
    print(h6.GetEntries() / n_tracks * 100, " - events with BS photon signal")
    h7 = gROOT.FindObject("h7")
    print(h7.GetEntries() / n_tracks * 100, " - events with BS e^+ signal")
    h8 = gROOT.FindObject("h8")
    print(h8.GetEntries() / n_tracks * 100, " - events with BS hadron signal")
    h9 = gROOT.FindObject("h9")
    print(h9.GetEntries() / n_tracks * 100, " - events with PS e^+ signal")
    h10 = gROOT.FindObject("h10")
    print(h10.GetEntries() / n_tracks * 100, " - events with PS hadron signal")
    h11 = gROOT.FindObject("h11")
    print(h11.GetEntries() / n_tracks * 100, " - events with BS any signal")
    h12 = gROOT.FindObject("h12")
    print(h12.GetEntries() / n_tracks * 100, " - events with PS any signal")

    input("wait")


def pass_2_trackers():
    # 0 - mixed
    # 1 - primary
    # 2- electrons
    # 3 - gammas
    # 4 - positrons
    # 5 - hadrons

    n_bs_hadr = 0
    n_bs_pos = 0
    n_bs_gamma = 0
    n_bs_e = 0
    n_ps_gamma = 0
    n_ps_e = 0
    n_ps_pos = 0
    n_ps_hadr = 0
    n_ps_prim = 0
    n_mix = 0

    n_both_bs_hadr = 0
    n_both_bs_pos = 0
    n_both_bs_gamma = 0
    n_both_bs_e = 0
    n_both_ps_gamma = 0
    n_both_ps_e = 0
    n_both_ps_pos = 0
    n_both_ps_hadr = 0
    n_both_ps_prim = 0
    n_both_mix = 0

    for idx, event in enumerate(t_mc_5gev):
        if idx % 10000 == 0:
            print(idx)
        for tr2_type, tr2_id, tr2_pz in zip(event.tr2_hit_type, event.tr2_hit_track_id, event.tr2_particle_pz):
            if (tr2_pz < 0) and (tr2_type == 5):
                n_bs_hadr += 1
                for tr1_type, tr1_id, tr1_pz in zip(event.tr1_hit_type, event.tr1_hit_track_id, event.tr1_particle_pz):
                    if (tr1_pz < 0) and (tr1_type == 5) and (tr1_id == tr2_id):
                        n_both_bs_hadr += 1

            if (tr2_pz < 0) and (tr2_type == 4):
                n_bs_pos += 1
                for tr1_type, tr1_id, tr1_pz in zip(event.tr1_hit_type, event.tr1_hit_track_id, event.tr1_particle_pz):
                    if (tr1_pz < 0) and (tr1_type == 4) and (tr1_id == tr2_id):
                        n_both_bs_pos += 1

            if (tr2_pz < 0) and (tr2_type == 3):
                n_bs_gamma += 1
                for tr1_type, tr1_id, tr1_pz in zip(event.tr1_hit_type, event.tr1_hit_track_id, event.tr1_particle_pz):
                    if (tr1_pz < 0) and (tr1_type == 3) and (tr1_id == tr2_id):
                        n_both_bs_gamma += 1

            if (tr2_pz < 0) and (tr2_type == 2):
                n_bs_e += 1
                for tr1_type, tr1_id, tr1_pz in zip(event.tr1_hit_type, event.tr1_hit_track_id, event.tr1_particle_pz):
                    if (tr1_pz < 0) and (tr1_type == 2) and (tr1_id == tr2_id):
                        n_both_bs_e += 1

            if (tr2_pz > 0) and (tr2_type == 3):
                n_ps_gamma += 1
                for tr1_type, tr1_id, tr1_pz in zip(event.tr1_hit_type, event.tr1_hit_track_id, event.tr1_particle_pz):
                    if (tr1_pz > 0) and (tr1_type == 3) and (tr1_id == tr2_id):
                        n_both_ps_gamma += 1

            if (tr2_pz > 0) and (tr2_type == 2):
                n_ps_e += 1
                for tr1_type, tr1_id, tr1_pz in zip(event.tr1_hit_type, event.tr1_hit_track_id, event.tr1_particle_pz):
                    if (tr1_pz > 0) and (tr1_type == 2) and (tr1_id == tr2_id):
                        n_both_ps_e += 1

            if (tr2_pz > 0) and (tr2_type == 1):
                n_ps_prim += 1
                for tr1_type, tr1_id, tr1_pz in zip(event.tr1_hit_type, event.tr1_hit_track_id, event.tr1_particle_pz):
                    if (tr1_pz > 0) and (tr1_type == 1) and (tr1_id == tr2_id):
                        n_both_ps_prim += 1

            if tr2_type == 0:
                n_mix += 1
                for tr1_type, tr1_id, tr1_pz in zip(event.tr1_hit_type, event.tr1_hit_track_id, event.tr1_particle_pz):
                    if tr1_type == 0 and (tr1_id == tr2_id):
                        n_both_mix += 1

            if (tr2_pz > 0) and (tr2_type == 5):
                n_ps_hadr += 1
                for tr1_type, tr1_id, tr1_pz in zip(event.tr1_hit_type, event.tr1_hit_track_id, event.tr1_particle_pz):
                    if (tr1_pz > 0) and (tr1_type == 5) and (tr1_id == tr2_id):
                        n_both_ps_hadr += 1

            if (tr2_pz > 0) and (tr2_type == 4):
                n_ps_pos += 1
                for tr1_type, tr1_id, tr1_pz in zip(event.tr1_hit_type, event.tr1_hit_track_id, event.tr1_particle_pz):
                    if (tr1_pz > 0) and (tr1_type == 4) and (tr1_id == tr2_id):
                        n_both_ps_pos += 1


    print(n_both_ps_prim / t_mc_5gev.GetEntries() * 100, " - events with primary particle")
    print(n_both_ps_e / t_mc_5gev.GetEntries() * 100, " - events with PS e^{-} signal")
    print(n_both_ps_gamma / t_mc_5gev.GetEntries() * 100, " - events with PS photon signal")
    print(n_both_mix / t_mc_5gev.GetEntries() * 100, " - events with mixed partcle signal")
    print(n_both_bs_e / t_mc_5gev.GetEntries() * 100, " - events with BS e^- signal")
    print(n_both_bs_gamma / t_mc_5gev.GetEntries() * 100, " - events with BS photon signal")
    print(n_both_bs_pos / t_mc_5gev.GetEntries() * 100, " - events with BS e^+ signal")
    print(n_both_bs_hadr / t_mc_5gev.GetEntries() * 100, " - events with BS hadron signal")
    print(n_both_ps_pos / t_mc_5gev.GetEntries() * 100, " - events with PS e^+ signal")
    print(n_both_ps_hadr / t_mc_5gev.GetEntries() * 100, " - events with PS hadron signal")

    print("RATIO numbers:")

    print(n_both_ps_prim / n_ps_prim * 100, " - events with primary particle", n_both_ps_prim, n_ps_prim)
    print(n_both_ps_e / n_ps_e * 100, " - events with PS e^{-} signal", n_both_ps_e, n_ps_e)
    print(n_both_ps_gamma / n_ps_gamma * 100, " - events with PS photon signal", n_both_ps_gamma, n_ps_gamma)
    print(n_both_mix / n_mix * 100, " - events with mixed partcle signal", n_both_mix, n_mix)
    print(n_both_bs_e / n_bs_e * 100, " - events with BS e^- signal", n_both_bs_e, n_bs_e)
    print(n_both_bs_gamma / n_bs_gamma * 100, " - events with BS photon signal", n_both_bs_gamma, n_bs_gamma)
    print(n_both_bs_pos / n_bs_pos * 100, " - events with BS e^+ signal", n_both_bs_pos, n_bs_pos)
    print(n_both_bs_hadr / n_bs_hadr * 100, " - events with BS hadron signal", n_both_bs_hadr, n_bs_hadr)
    print(n_both_ps_pos / n_ps_pos * 100, " - events with PS e^+ signal", n_both_ps_pos, n_ps_pos)
    print(n_both_ps_hadr / n_ps_hadr * 100, " - events with PS hadron signal", n_both_ps_hadr, n_ps_hadr)

    input("wait")


def beam_shape():
    t_data_5gev.Draw("cal_cluster_y[0]>>h_data(64, 80., 195.2)")
    t_mc_5gev.Draw("cal_cluster_y[0]>>h_mc(64, 80., 195.2)")

    h_data = gROOT.FindObject("h_data")
    h_data.Scale(1./t_data_5gev.GetEntries())
    h_mc = gROOT.FindObject("h_mc")
    h_mc.Scale(1./t_mc_5gev.GetEntries())
    h_data.Divide(h_mc)

    h_data.Draw()
    for i in range(64):
        print("Bin {}".format(i), h_data.GetBinContent(i))
    # h_mc.Draw("histosame")
    input('wait')


def check_clustering():
    # Position of the cluster and its energy should correlate
    # t_data_5gev.Draw("cal_cluster_y[0]:cal_cluster_energy[0]>>h1(100,0.,500,64,80.,195.2)", "", "colz")
    # input("meow")
    # Energy of the e+gamma clusters should correlate and be a constant summ.
    t_data_5gev.Draw("cal_cluster_energy[1]:cal_cluster_energy[0]>>h2(500,0, 500, 250, 0, 250)", "cal_n_clusters > 1", "colz")
    input("woof")
    # Position of the cluster if it is reconstructed as 1
    t_data_5gev.Draw("cal_cluster_y:cal_n_clusters>>h2(3,0,3, 88, 116., 195.2)", "", "colz")
    input("quack")


def check_event():
    for i, ev in enumerate(t_data_5gev):
        if ev.cal_n_clusters != 1:
            continue
        print("Clusters:", ev.cal_n_clusters, "in event:", i)
        t_data_5gev.Draw("cal_hit_pad:cal_hit_sector>>h1(4,0.,4,64,0,64)", "cal_hit_energy*(Entry$ == {})".format(i), "colz")
        input("check")


def identification():
    cut_off = np.arange(0.5, 30.1, 0.1)

    n_gen = 0
    n_gen_reco_e = np.zeros_like(cut_off)
    n_reco_e = np.zeros_like(cut_off)
    n_gen_reco_ph = np.zeros_like(cut_off)
    n_reco_ph = np.zeros_like(cut_off)


    n_events = t_data_5gev.GetEntries()
    for i, event in enumerate(t_data_5gev):
        if i % 1000 == 0:
            print(i, "event out of", n_events)
            if i == 20000:
                break
        # Event with 2 clusters. 1 electron 1 photon.
        if not (event.cal_n_clusters == 2
                and 130 < event.cal_cluster_y[0] < 150
                and 158 < event.cal_cluster_y[1] < 178
                and event.tr1_n_hits > 0
                and event.tr2_n_hits > 0):
            continue
        n_gen += 1

        cal_e_y = event.cal_cluster_y[0]
        cal_gamma_y = event.cal_cluster_y[1]
        tr1y = np.array(event.tr1_hit_y)
        tr2y = np.array(event.tr2_hit_y)
        # print("reco_e", n_reco_e)
        # print("reco_ph", n_reco_e)
        for idx, cut in enumerate(cut_off):
            # If electron cluster has hits in both trackers within cut. Therefore it is reco as e
            if any(abs(cal_e_y - tr1y) < cut) and any(abs(cal_e_y - tr2y) < cut):
                n_gen_reco_e[idx] += 1
                n_reco_e[idx] += 1
            # If photon cluster has hits in both trackers within cut. Therefore it is reco as e
            if any(abs(cal_gamma_y - tr1y) < cut) and any(abs(cal_gamma_y - tr2y) < cut):
                n_reco_e[idx] += 1

            # If electron cluster has no hits in both trackers within cut. Therefore it is reco as photon
            if (not any(abs(cal_e_y - tr1y) < cut)) and (not any(abs(cal_e_y - tr2y) < cut)):
                n_reco_ph[idx] += 1

            # If gamma cluster has no hits in both trackers within cut. Therefore it is reco as photon
            if (not any(abs(cal_gamma_y - tr1y) < cut)) and (not any(abs(cal_gamma_y - tr2y) < cut)):
                n_gen_reco_ph[idx] += 1
                n_reco_ph[idx] += 1

    # Get the graph
    eff_e = n_gen_reco_e / n_gen * 100
    purity_e = n_gen_reco_e / n_reco_e * 100

    eff_ph = n_gen_reco_ph / n_gen * 100
    purity_ph = n_gen_reco_ph / n_reco_ph * 100

    gr_e = TGraph(len(cut_off), eff_e, purity_e)
    gr_ph = TGraph(len(cut_off), eff_ph, purity_ph)
    gr_ph.SetLineColor(2)
    gr_ph.SetMarkerColor(2)
    gr_ph.SetMarkerStyle(22)

    gr_e.GetXaxis().SetTitle("Efficiency, %")
    gr_e.GetYaxis().SetTitle("Purity, %")
    gr_e.SetTitle("Electron ident algo")
    gr_ph.SetTitle("Photon ident algo")

    gr_e.Draw("APL")
    gr_ph.Draw("samePL")
    input("wait")
    # gr.Write("puref")














# check_event()
# check_clustering()
# identification()
# y_distance()
# dep_energy()
# particle_energy()
# track_len()
calculate_numbers()
