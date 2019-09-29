from ROOT import TFile, gROOT, TGraphErrors, TGraph, TF1, gStyle, TF1, nullptr, TH2F, TMath
from ROOT import TCanvas, TColor, TGaxis, TH1F, TPad
from ROOT import kBlack, kBlue, kRed

import numpy as np

gROOT.SetStyle('ATLAS')
gROOT.SetBatch(0)

# Data trees
f_data_5gev = TFile.Open("../extracted/extracted_data_5gev.root", 'read')
f_data_4gev = TFile.Open('../extracted/extracted_data_4gev.root', 'read')
f_data_3gev = TFile.Open('../extracted/extracted_data_3gev.root', 'read')
f_data_2gev = TFile.Open('../extracted/extracted_data_2gev.root', 'read')
f_data_1gev = TFile.Open('../extracted/extracted_data_1gev.root', 'read')
t_data_5gev = f_data_5gev.data
t_data_4gev = f_data_4gev.data
t_data_3gev = f_data_3gev.data
t_data_2gev = f_data_2gev.data
t_data_1gev = f_data_1gev.data

# MC trees
f_mc_5gev = TFile.Open("../extracted/extracted_lucas_5gev.root", 'read')
# f_mc_4gev = TFile.Open("../extracted/extracted_lucas_4gev.root", 'read')
# f_mc_3gev = TFile.Open("../extracted/extracted_lucas_3gev.root", 'read')
# f_mc_2gev = TFile.Open("../extracted/extracted_lucas_2gev.root", 'read')
# f_mc_1gev = TFile.Open("../extracted/extracted_lucas_1gev.root", 'read')
t_mc_5gev = f_mc_5gev.mc
# t_mc_4gev = f_mc_4gev.mc
# t_mc_3gev = f_mc_3gev.mc
# t_mc_2gev = f_mc_2gev.mc
# t_mc_1gev = f_mc_1gev.mc


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


# Part for fancy ratio plots
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
# Part for fancy ratio plots


def method1():
    """Fit the horizontal line through reconstructed shower and hit in tracker2.
    Calculate expected position of hit in tracker1 from the fit of those 2 points.
    Plot histo of residuals of expected and measured positions.
    Plot efficiency and purity points for various cut-off values on the residual."""
    h_reco_gen = TH1F("h_reco_gen", "title", 800, -20., 20.)
    h_reco = TH1F("h_reco", "title", 800, -20., 20.)
    n_gen = 0

    n_events = t_mc_5gev.GetEntries()
    x = np.array([6., 24.])
    y_err = np.array([0.52, 0.44])

    for i, event in enumerate(t_mc_5gev):
        if i % 1000 == 0:
            print(i, "event out of", n_events)
        if event.cal_n_clusters == 0:
            continue

        for y1, p1 in zip(event.tr1_hit_y, event.tr1_hit_is_prime):
            for y2, p2 in zip(event.tr2_hit_y, event.tr2_hit_is_prime):

                # This is algorithm
                y = np.array([y2, event.cal_cluster_y[0]])
                gr_track = TGraphErrors(2, x, y, nullptr, y_err)
                gr_track.Fit("pol0", "Q")
                residual = y1 - gr_track.GetFunction("pol0").Eval(1.)

                # there was a track
                if p1 == 1 and p2 == 1:
                    n_gen += 1
                    h_reco_gen.Fill(residual - 0.225)

                h_reco.Fill(residual - 0.225)

    # Get the graph
    cut_off = np.arange(0., 20., 0.1)
    eff = []
    purity = []

    for cut in cut_off:
        bmin = h_reco_gen.GetXaxis().FindBin(-cut)
        bmax = h_reco_gen.GetXaxis().FindBin(cut)
        n_reco_gen = h_reco_gen.Integral(bmin, bmax)
        n_reco = h_reco.Integral(bmin, bmax)

        if n_gen != 0 and n_reco != 0:
            eff.append(n_reco_gen / n_gen * 100)
            purity.append(n_reco_gen / n_reco * 100)

    gr = TGraph(len(cut_off), np.array(eff), np.array(purity))
    gr.GetXaxis().SetTitle("Efficiency, %")
    gr.GetYaxis().SetTitle("Purity, %")
    gr.SetTitle("Tr2 Cal fit")

    return gr


def method2():
    """Fit the horizontal line through hits in tracker1 and tracker2.
    Calculate expected position of hit in the calorimeter from the fit of those 2 points.
    Plot histo of residuals of expected and measured positions of those hits.
    Plot efficiency and purity points for various cut-off values on the residual."""
    h_reco_gen = TH1F("h_reco_gen", "title", 800, -20., 20.)
    h_reco = TH1F("h_reco", "title", 800, -20., 20.)
    n_gen = 0

    n_events = t_mc_5gev.GetEntries()
    x = np.array([1., 6.])
    y_err = np.array([0.52, 0.52])

    for i, event in enumerate(t_mc_5gev):
        if i % 1000 == 0:
            print(i, "event out of", n_events)
        if event.cal_n_clusters == 0:
            continue

        for y1, p1 in zip(event.tr1_hit_y, event.tr1_hit_is_prime):
            for y2, p2 in zip(event.tr2_hit_y, event.tr2_hit_is_prime):

                # This is algorithm
                y = np.array([y1, y2])
                gr_track = TGraphErrors(2, x, y, nullptr, y_err)
                gr_track.Fit("pol0", "Q")
                residual = event.cal_cluster_y[0] - gr_track.GetFunction("pol0").Eval(24.)

                # there was a track
                if p1 == 1 and p2 == 1:
                    n_gen += 1
                    h_reco_gen.Fill(residual)

                h_reco.Fill(residual)

    # Get the graph
    cut_off = np.arange(0., 20., 0.1)
    eff = []
    purity = []

    for cut in cut_off:
        bmin = h_reco_gen.GetXaxis().FindBin(-cut)
        bmax = h_reco_gen.GetXaxis().FindBin(cut)
        n_reco_gen = h_reco_gen.Integral(bmin, bmax)
        n_reco = h_reco.Integral(bmin, bmax)

        if n_gen != 0 and n_reco != 0:
            eff.append(n_reco_gen / n_gen * 100)
            purity.append(n_reco_gen / n_reco * 100)

    gr = TGraph(len(cut_off), np.array(eff), np.array(purity))
    gr.GetXaxis().SetTitle("Efficiency, %")
    gr.GetYaxis().SetTitle("Purity, %")
    gr.SetTitle("Tr1 Tr2 fit")

    return gr


def method3():
    """Fit the horizontal line through hits in tracker1, tracker2 and calorimeter.
    Calculate expected position of hit in the calorimeter from the fit of those 3 points.
    Plot histo of residuals of expected and measured positions of the shower.
    Plot efficiency and purity points for various cut-off values on the residual."""
    h_reco_gen = TH1F("h_reco_gen", "title", 800, -20., 20.)
    h_reco = TH1F("h_reco", "title", 800, -20., 20.)
    n_gen = 0

    n_events = t_mc_5gev.GetEntries()
    x = np.array([1., 6., 24])
    y_err = np.array([0.52, 0.52, 0.44])

    for i, event in enumerate(t_mc_5gev):
        if i % 1000 == 0:
            print(i, "event out of", n_events)
        if event.cal_n_clusters == 0:
            continue

        for y1, p1 in zip(event.tr1_hit_y, event.tr1_hit_is_prime):
            for y2, p2 in zip(event.tr2_hit_y, event.tr2_hit_is_prime):

                # This is algorithm
                y = np.array([y1, y2, event.cal_cluster_y[0]])
                gr_track = TGraphErrors(3, x, y, nullptr, y_err)
                gr_track.Fit("pol0", "Q")
                residual = event.cal_cluster_y[0] - gr_track.GetFunction("pol0").Eval(24.)

                # there was a track
                if p1 == 1 and p2 == 1:
                    n_gen += 1
                    h_reco_gen.Fill(residual)

                h_reco.Fill(residual)

    # Get the graph
    cut_off = np.arange(0., 20., 0.1)
    eff = []
    purity = []

    for cut in cut_off:
        bmin = h_reco_gen.GetXaxis().FindBin(-cut)
        bmax = h_reco_gen.GetXaxis().FindBin(cut)
        n_reco_gen = h_reco_gen.Integral(bmin, bmax)
        n_reco = h_reco.Integral(bmin, bmax)

        if n_gen != 0 and n_reco != 0:
            eff.append(n_reco_gen / n_gen * 100)
            purity.append(n_reco_gen / n_reco * 100)

    gr = TGraph(len(cut_off), np.array(eff), np.array(purity))
    gr.GetXaxis().SetTitle("Efficiency, %")
    gr.GetYaxis().SetTitle("Purity, %")
    gr.SetTitle("Tr1 Tr2 Cal fit")

    return gr



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


gr1 = method1()
gr1.Draw("AP")

gr2 = method2()
gr2.SetMarkerColor(2)
gr2.SetMarkerStyle(21)
gr2.Draw("Psame")

gr3 = method3()
gr3.SetMarkerColor(4)
gr3.SetMarkerStyle(23)
gr3.Draw("Psame")

input("wait")
# pad_spectra()
