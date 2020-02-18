from ROOT import TFile, gROOT, TCanvas, TPad, kBlack

gROOT.SetStyle('ATLAS')
# Set arg to 0 to print pictures fast
gROOT.SetBatch(0)

# Files
f_data_5gev_nomerge = TFile.Open("./extracted_trees/extracted_data_5gev_nomerge.root", 'read')
f_data_5gev = TFile.Open("./extracted_trees/extracted_data_5gev.root", 'read')
f_data_4gev = TFile.Open('./extracted_trees/extracted_data_4gev.root', 'read')
f_data_3gev = TFile.Open('./extracted_trees/extracted_data_3gev.root', 'read')
f_data_2gev = TFile.Open('./extracted_trees/extracted_data_2gev.root', 'read')
f_data_1gev = TFile.Open('./extracted_trees/extracted_data_1gev.root', 'read')
f_mc_5gev = TFile.Open("./extracted_trees/extracted_mc_5gev.root", 'read')
# f_mc_4gev = TFile.Open("./extracted_trees/extracted_lucas_4gev.root", 'read')
# f_mc_3gev = TFile.Open("./extracted_trees/extracted_lucas_3gev.root", 'read')
# f_mc_2gev = TFile.Open("./extracted_trees/extracted_lucas_2gev.root", 'read')
# f_mc_1gev = TFile.Open("./extracted_trees/extracted_lucas_1gev.root", 'read')

# Trees
t_data_5gev_nomerge = f_data_5gev_nomerge.data
t_data_5gev = f_data_5gev.data
t_data_4gev = f_data_4gev.data
t_data_3gev = f_data_3gev.data
t_data_2gev = f_data_2gev.data
t_data_1gev = f_data_1gev.data
t_mc_5gev = f_mc_5gev.mc
# t_mc_4gev = f_mc_4gev.mc
# t_mc_3gev = f_mc_3gev.mc
# t_mc_2gev = f_mc_2gev.mc
# t_mc_1gev = f_mc_1gev.mc


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


def plot_n_clusters(opt_print=0):
    t_data_5gev.Draw("cal_n_clusters>>h_data(6,0,6)")
    h_data = gROOT.FindObject("h_data")
    h_data.SetTitle("Data")
    h_data.GetXaxis().SetTitle("N_{cal,clusters}")
    h_data.GetYaxis().SetTitle("N_{clusters}/N_{events}")
    h_data.SetLineWidth(3)
    h_data.Scale(1. / t_data_5gev.GetEntries())

    t_mc_5gev.Draw("cal_n_clusters>>h_mc(6,0,6)")
    h_mc = gROOT.FindObject("h_mc")
    h_mc.SetTitle("MC")
    h_mc.SetLineColor(2)
    h_mc.SetLineWidth(3)
    h_mc.Scale(1. / t_mc_5gev.GetEntries())

    c, pad1, pad2 = createCanvasPads("plot_n_clusters")
    pad1.cd()
    h_data.Draw("histo")
    h_mc.Draw("histo same")
    h_data.GetYaxis().SetRangeUser(0., 1.)
    pad1.BuildLegend()

    h_ratio = createRatio(h_data, h_mc)
    pad2.cd()
    h_ratio.Draw("EP")
    h_ratio.GetYaxis().UnZoom()

    if opt_print == 0:
        input("wait")
    elif opt_print == 1:
        c.Print("./plot_n_clusters.png")


def plot_energy(opt_print=0):
    t_data_5gev.Draw("cal_cluster_energy*0.0885>>h_data(500,0.,60)")
    h_data = gROOT.FindObject("h_data")
    h_data.SetTitle("Data")
    h_data.GetXaxis().SetTitle("E_{cal,cluster}, MeV")
    h_data.GetYaxis().SetTitle("N_{clusters}/N_{events}")
    h_data.SetLineWidth(3)
    h_data.Scale(1. / t_data_5gev.GetEntries())

    t_mc_5gev.Draw("cal_cluster_energy*0.0885>>h_mc(500,0.,60)")
    h_mc = gROOT.FindObject("h_mc")
    h_mc.SetTitle("MC")
    h_mc.SetLineColor(2)
    h_mc.SetLineWidth(3)
    h_mc.Scale(1. / t_mc_5gev.GetEntries())

    c, pad1, pad2 = createCanvasPads("plot_energy")
    pad1.cd()
    h_data.Draw("histo")
    h_mc.Draw("histo same")
    h_data.GetYaxis().SetRangeUser(0., 1.05 * max(h_data.GetMaximum(), h_mc.GetMaximum()))
    pad1.BuildLegend()

    h_ratio = createRatio(h_data, h_mc)
    pad2.cd()
    h_ratio.Draw("EP")
    h_ratio.GetYaxis().UnZoom()

    if opt_print == 0:
        input("wait")
    elif opt_print == 1:
        c.Print("./plot_energy.png")


def plot_y(opt_print=0):
    t_data_5gev.Draw("cal_cluster_y>>h_data(44,116.,195.2)")
    h_data = gROOT.FindObject("h_data")
    h_data.SetTitle("Data")
    h_data.GetXaxis().SetTitle("y_{cal,cluster}, mm")
    h_data.GetYaxis().SetTitle("N_{clusters}/N_{events}")
    h_data.SetLineWidth(3)
    h_data.Scale(1. / t_data_5gev.GetEntries())

    t_mc_5gev.Draw("cal_cluster_y>>h_mc(44,116.,195.2)")
    h_mc = gROOT.FindObject("h_mc")
    h_mc.SetTitle("MC")
    h_mc.SetLineColor(2)
    h_mc.SetLineWidth(3)
    h_mc.Scale(1. / t_mc_5gev.GetEntries())

    c, pad1, pad2 = createCanvasPads("plot_y")
    pad1.cd()
    h_data.Draw("histo")
    h_mc.Draw("histo same")
    h_data.GetYaxis().SetRangeUser(0., 1.05 * max(h_data.GetMaximum(), h_mc.GetMaximum()))
    pad1.BuildLegend()

    h_ratio = createRatio(h_data, h_mc)
    pad2.cd()
    h_ratio.Draw("EP")
    h_ratio.GetYaxis().UnZoom()

    if opt_print == 0:
        input("wait")
    elif opt_print == 1:
        c.Print("./plot_y.png")


def plot_x(opt_print=0):
    t_data_5gev.Draw("cal_cluster_x>>h_data(100, -15., 15.)")
    h_data = gROOT.FindObject("h_data")
    h_data.SetTitle("Data")
    h_data.GetXaxis().SetTitle("x_{cal,cluster}, mm")
    h_data.GetYaxis().SetTitle("N_{clusters}/N_{events}")
    h_data.SetLineWidth(3)
    h_data.Scale(1. / t_data_5gev.GetEntries())

    t_mc_5gev.Draw("cal_cluster_x>>h_mc(100, -15., 15.)")
    h_mc = gROOT.FindObject("h_mc")
    h_mc.SetTitle("MC")
    h_mc.SetLineColor(2)
    h_mc.SetLineWidth(3)
    h_mc.Scale(1. / t_mc_5gev.GetEntries())

    c, pad1, pad2 = createCanvasPads("plot_x")
    pad1.cd()
    h_data.Draw("histo")
    h_mc.Draw("histo same")
    h_data.GetYaxis().SetRangeUser(0., 1.05 * max(h_data.GetMaximum(), h_mc.GetMaximum()))
    pad1.BuildLegend()

    h_ratio = createRatio(h_data, h_mc)
    pad2.cd()
    h_ratio.Draw("EP")
    h_ratio.GetYaxis().UnZoom()

    if opt_print == 0:
        input("wait")
    elif opt_print == 1:
        c.Print("./plot_x.png")


# plot_n_clusters()
plot_energy()
# plot_x()
# plot_y()
