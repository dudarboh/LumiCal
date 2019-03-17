from ROOT import TFile, gROOT, TGraphErrors, TH1F, TGraph, TH2F, TCanvas, TPad, TF1, gStyle, TColor
import numpy as np
import array


class Data():
    def __init__(self):
        self.file_data = TFile.Open('./trees/extracted_data_merge_on.root')
        self.tree_data = self.file_data.data
        self.tree_data.AddFriend('nomerge = data', './trees/extracted_data_merge_off.root')

        self.file_mc = TFile.Open('./trees/extracted_mc_merge_on.root')
        self.tree_mc = self.file_mc.mc
        self.tree_mc.AddFriend('nomerge = mc', './trees/extracted_mc_merge_off.root')

        self.output_file = TFile('RENAME.root', 'recreate')
        self.output_file.cd()
        self.n_events_data = self.tree_data.GetEntries()
        self.n_events_mc = self.tree_mc.GetEntries()

    # Newly made for presentation
    def calorimeter_plot_y(self):
        canvas = TCanvas('calorimeter_plot_y', 'title', 1024, 768)
        n_bins = 180
        first = 110
        last = 200

        self.tree_data.Draw('cal_cluster_y[0]>>h_data({}, {}, {})'.format(n_bins, first, last))
        self.tree_mc.Draw('cal_cluster_y[0]>>h_mc({}, {}, {})'.format(n_bins, first, last))
        self.tree_mc.Draw('cal_true_hit_y>>h_true({}, {}, {})'.format(n_bins, first, last))
        h_data = gROOT.FindObject('h_data')
        h_mc = gROOT.FindObject('h_mc')
        h_true = gROOT.FindObject('h_true')

        h_mc.Draw('histo')
        h_mc.Scale(1. / h_mc.GetEntries())
        h_data.Draw('histosame')
        h_data.Scale(1. / h_data.GetEntries())
        h_true.Draw('histosame')
        h_true.Scale(1. / h_mc.GetEntries())

        h_mc.SetMaximum(max([h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin()), h_true.GetBinContent(h_true.GetMaximumBin())]))

        h_mc.SetFillColor(5)
        h_mc.SetTitle('MC;y, [mm];#frac{N_{ev}}{N_{tot}}')

        h_data.SetLineWidth(3)
        h_data.SetTitle('Data')

        h_true.SetLineWidth(3)
        h_true.SetLineColor(2)
        h_true.SetTitle('Generated')

        canvas.BuildLegend()
        h_mc.SetTitle('Shower y position')
        gStyle.SetOptStat(1110)
        canvas.Write('calorimeter_plot_y')

    def calorimeter_plot_x(self):
        canvas = TCanvas('calorimeter_plot_x', 'title', 1024, 768)
        n_bins = 60
        first = -15
        last = 15

        self.tree_data.Draw('cal_cluster_x[0]>>h_data({}, {}, {})'.format(n_bins, first, last))
        self.tree_mc.Draw('cal_cluster_x[0]>>h_mc({}, {}, {})'.format(n_bins, first, last))
        self.tree_mc.Draw('cal_true_hit_x>>h_true({}, {}, {})'.format(n_bins, first, last))
        h_data = gROOT.FindObject('h_data')
        h_mc = gROOT.FindObject('h_mc')
        h_true = gROOT.FindObject('h_true')

        h_mc.Draw('histo')
        h_mc.Scale(1. / h_mc.GetEntries())
        h_data.Draw('histosame')
        h_data.Scale(1. / h_data.GetEntries())
        h_true.Draw('histosame')
        h_true.Scale(1. / h_mc.GetEntries())

        h_mc.SetMaximum(max([h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin()), h_true.GetBinContent(h_true.GetMaximumBin())]))

        h_mc.SetFillColor(5)
        h_mc.SetTitle('MC;x, [mm];#frac{N_{ev}}{N_{tot}}')

        h_data.SetLineWidth(3)
        h_data.SetTitle('Data')

        h_true.SetLineWidth(3)
        h_true.SetLineColor(2)
        h_true.SetTitle('Generated')

        canvas.BuildLegend()
        h_mc.SetTitle('Shower x position')
        gStyle.SetOptStat(1110)
        canvas.Write('calorimeter_plot_x')

    def calorimeter_plot_energy(self):
        canvas = TCanvas('calorimeter_plot_energy', 'title', 1024, 768)
        n_bins = 900
        first = 0
        last = 900

        self.tree_data.Draw('cal_cluster_energy[0]>>h_data({}, {}, {})'.format(n_bins, first, last))
        self.tree_mc.Draw('cal_cluster_energy[0]>>h_mc({}, {}, {})'.format(n_bins, first, last))
        h_data = gROOT.FindObject('h_data')
        h_mc = gROOT.FindObject('h_mc')

        h_mc.Draw('histo')
        h_mc.Scale(1. / h_mc.GetEntries())
        h_data.Draw('histosame')
        h_data.Scale(1. / h_data.GetEntries())

        h_mc.SetMaximum(max(h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin())))

        h_mc.SetFillColor(5)
        h_mc.SetTitle('MC;energy, [MIP];#frac{N_{ev}}{N_{tot}}')

        h_data.SetLineWidth(3)
        h_data.SetTitle('Data')

        canvas.BuildLegend()
        h_mc.SetTitle('Shower energy')
        gStyle.SetOptStat(1110)
        canvas.Write('calorimeter_plot_energy')

    def calorimeter_plot_n_pads(self):
        canvas = TCanvas('calorimeter_plot_n_pads', 'title', 1024, 768)
        n_bins = 80
        first = 0
        last = 80

        self.tree_data.Draw('cal_cluster_n_pads[0]>>h_data({}, {}, {})'.format(n_bins, first, last))
        self.tree_mc.Draw('cal_cluster_n_pads[0]>>h_mc({}, {}, {})'.format(n_bins, first, last))
        h_data = gROOT.FindObject('h_data')
        h_mc = gROOT.FindObject('h_mc')

        h_mc.Draw('histo')
        h_mc.Scale(1. / h_mc.GetEntries())
        h_data.Draw('histosame')
        h_data.Scale(1. / h_data.GetEntries())

        h_mc.SetMaximum(max(h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin())))

        h_mc.SetFillColor(5)
        h_mc.SetTitle('MC;N_{pads};#frac{N_{ev}}{N_{tot}}')

        h_data.SetLineWidth(3)
        h_data.SetTitle('Data')

        canvas.BuildLegend()
        h_mc.SetTitle('Number of hits in shower')
        gStyle.SetOptStat(1110)
        canvas.Write('calorimeter_plot_n_pads')

    def calorimeter_plot_n_clusters(self):
        canvas = TCanvas('calorimeter_plot_n_clusters', 'title', 1024, 768)
        n_bins = 7
        first = 0
        last = 7

        self.tree_data.Draw('cal_n_clusters>>h_data({}, {}, {})'.format(n_bins, first, last))
        self.tree_mc.Draw('cal_n_clusters>>h_mc({}, {}, {})'.format(n_bins, first, last))
        h_data = gROOT.FindObject('h_data')
        h_mc = gROOT.FindObject('h_mc')

        h_mc.Draw('histo')
        h_mc.Scale(1. / h_mc.GetEntries())
        h_data.Draw('histosame')
        h_data.Scale(1. / h_data.GetEntries())

        h_mc.SetMaximum(max(h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin())))

        h_mc.SetFillColor(5)
        h_mc.SetTitle('MC;N_{clusters};#frac{N_{ev}}{N_{tot}}')

        h_data.SetLineWidth(3)
        h_data.SetTitle('Data')

        canvas.BuildLegend()
        h_mc.SetTitle('Number of clusters')
        gStyle.SetOptStat(1110)
        canvas.Write('calorimeter_plot_n_clusters in event')

    def calorimeter_plot_layer(self):
        canvas = TCanvas('calorimeter_plot_layer', 'title', 1024, 768)
        n_bins = 70
        first = 0
        last = 7

        self.tree_data.Draw('cal_cluster_layer[0]>>h_data({}, {}, {})'.format(n_bins, first, last))
        self.tree_mc.Draw('cal_cluster_layer[0]>>h_mc({}, {}, {})'.format(n_bins, first, last))
        h_data = gROOT.FindObject('h_data')
        h_mc = gROOT.FindObject('h_mc')

        h_mc.Draw('histo')
        h_mc.Scale(1. / h_mc.GetEntries())
        h_data.Draw('histosame')
        h_data.Scale(1. / h_data.GetEntries())

        h_mc.SetMaximum(max(h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin())))

        h_mc.SetFillColor(5)
        h_mc.SetTitle('MC;layer;#frac{N_{ev}}{N_{tot}}')

        h_data.SetLineWidth(3)
        h_data.SetTitle('Data')

        canvas.BuildLegend()
        h_mc.SetTitle('Shower layer')
        gStyle.SetOptStat(1110)
        canvas.Write('calorimeter_plot_layer')

    def calorimeter_plot_n_pads_vs_energy(self):
        n_bins1 = 200
        first1 = 0
        last1 = 600
        n_bins2 = 90
        first2 = 0
        last2 = 90

        self.tree_data.Draw('cal_cluster_n_pads[0]:cal_cluster_energy[0]>>h_data({}, {}, {}, {}, {}, {})'.format(n_bins1, first1, last1, n_bins2, first2, last2))
        h_data = gROOT.FindObject('h_data')
        h_data.Write('calorimeter_plot_n_pads_vs_energy')

        self.tree_mc.Draw('cal_cluster_n_pads[0]:cal_cluster_energy[0]>>h_mc({}, {}, {}, {}, {}, {})'.format(n_bins1, first1, last1, n_bins2, first2, last2))
        h_mc = gROOT.FindObject('h_mc')
        h_mc.Write('calorimeter_plot_n_pads_vs_energy_mc')

    def tr1_distance_to_shower(self):
        canvas = TCanvas('tr1_distance_to_shower', 'title', 1024, 768)
        n_bins = 200
        first = -6
        last = 6

        self.tree_data.Draw('tr1_cluster_y[0]-cal_cluster_y[0]>>h_data({}, {}, {})'.format(n_bins, first, last))
        self.tree_mc.Draw('tr1_cluster_y[0]-cal_cluster_y[0]>>h_mc({}, {}, {})'.format(n_bins, first, last))
        h_data = gROOT.FindObject('h_data')
        h_mc = gROOT.FindObject('h_mc')

        h_mc.Draw('histo')
        h_mc.Scale(1. / h_mc.GetEntries())
        h_data.Draw('histosame')
        h_data.Scale(1. / h_data.GetEntries())

        h_mc.SetMaximum(max(h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin())))

        h_mc.SetFillColor(5)
        h_mc.SetTitle('MC;d_{cluster-shower}, [mm];#frac{N_{ev}}{N_{tot}}')

        h_data.SetLineWidth(3)
        h_data.SetTitle('Data')

        canvas.BuildLegend()
        h_mc.SetTitle('Distance:Tracker1 main cluster to the shower')
        gStyle.SetOptStat(1110)
        canvas.Write('tr1_distance_to_shower')

    def tr2_distance_to_shower(self):
        canvas = TCanvas('tr2_distance_to_shower', 'title', 1024, 768)
        n_bins = 200
        first = -6
        last = 6

        self.tree_data.Draw('tr2_cluster_y[0]-cal_cluster_y[0]>>h_data({}, {}, {})'.format(n_bins, first, last))
        self.tree_mc.Draw('tr2_cluster_y[0]-cal_cluster_y[0]>>h_mc({}, {}, {})'.format(n_bins, first, last))
        h_data = gROOT.FindObject('h_data')
        h_mc = gROOT.FindObject('h_mc')

        h_mc.Draw('histo')
        h_mc.Scale(1. / h_mc.GetEntries())
        h_data.Draw('histosame')
        h_data.Scale(1. / h_data.GetEntries())

        h_mc.SetMaximum(max(h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin())))

        h_mc.SetFillColor(5)
        h_mc.SetTitle('MC;d_{cluster-shower}, [mm];#frac{N_{ev}}{N_{tot}}')

        h_data.SetLineWidth(3)
        h_data.SetTitle('Data')

        canvas.BuildLegend()
        h_mc.SetTitle('Distance:Tracker2 main cluster to the shower')
        gStyle.SetOptStat(1110)
        canvas.Write('tr2_distance_to_shower')

    def tr1_efficiency(self):
        canvas = TCanvas('tr1_efficiency', 'title', 1024, 768)

        general_cuts = 'cal_n_clusters == 1 && tr2_n_clusters == 1 && cal_cluster_y[0] > 158 && cal_cluster_y[0] < 168\
                && abs(tr2_cluster_y - cal_cluster_y[0]) < 1.8'

        self.tree_data.Draw('>>h_data_ev_list', general_cuts)
        self.tree_mc.Draw('>>h_mc_ev_list', general_cuts)
        h_data_ev_list = gROOT.FindObject('h_data_ev_list')
        h_mc_ev_list = gROOT.FindObject('h_mc_ev_list')

        total_entries_data = h_data_ev_list.GetN()
        total_entries_mc = h_mc_ev_list.GetN()

        x_arr = np.arange(0, 5, 0.05)

        efficiency_data = np.zeros_like(x_arr)
        efficiency_mc = np.zeros_like(x_arr)

        for i, x in enumerate(x_arr):
            selection_cuts = 'cal_n_clusters == 1 && tr2_n_clusters == 1 && cal_cluster_y[0] > 158 && cal_cluster_y[0] < 168\
                    && abs(tr2_cluster_y - cal_cluster_y[0]) < 1.8 && abs(tr1_cluster_y - cal_cluster_y[0]) < {}'.format(x)

            self.tree_data.Draw('>>h_data_ev_list', selection_cuts)
            self.tree_mc.Draw('>>h_mc_ev_list', selection_cuts)
            h_data_ev_list = gROOT.FindObject('h_data_ev_list')
            h_mc_ev_list = gROOT.FindObject('h_mc_ev_list')

            efficiency_data[i] = h_data_ev_list.GetN() / total_entries_data * 100
            efficiency_mc[i] = h_mc_ev_list.GetN() / total_entries_mc * 100

        efficiency_data_graph = TGraph(len(x_arr), x_arr, efficiency_data)
        efficiency_data_graph.SetTitle('Data')

        efficiency_data_graph.Draw()

        efficiency_mc_graph = TGraph(len(x_arr), x_arr, efficiency_mc)
        efficiency_mc_graph.SetTitle('MC')
        efficiency_mc_graph.SetLineColor(2)
        efficiency_mc_graph.SetMarkerColor(2)
        efficiency_mc_graph.Draw('samePL')

        canvas.BuildLegend()
        eff_data_graph_copy = efficiency_data_graph.Clone()
        eff_mc_graph_copy = efficiency_mc_graph.Clone()
        efficiency_data_graph.SetTitle('Tracker1 efficiency;|y_{tr1} - y_{shower}|;Efficiency, %')

        # Scaled subpad

        subpad = TPad("subpad", "", 0.35, 0.35, 0.9, 0.8)
        subpad.Draw()
        subpad.cd()
        subpad.SetGridx()
        subpad.SetGridy()

        eff_data_graph_copy.Draw()
        eff_data_graph_copy.SetMarkerStyle(7)

        eff_data_graph_copy.SetMaximum(100.5)
        eff_data_graph_copy.SetMinimum(98)
        eff_data_graph_copy.GetXaxis().SetLimits(1, 5)
        eff_data_graph_copy.SetTitle('')

        eff_mc_graph_copy.Draw('sameLP')
        eff_mc_graph_copy.SetMarkerStyle(7)

        canvas.Write('tr1_efficiency')

    def tr2_efficiency(self):
        canvas = TCanvas('tr2_efficiency', 'title', 1024, 768)

        general_cuts = 'cal_n_clusters == 1 && tr1_n_clusters == 1 && cal_cluster_y[0] > 158 && cal_cluster_y[0] < 168\
                && abs(tr1_cluster_y - cal_cluster_y[0]) < 1.8'

        self.tree_data.Draw('>>h_data_ev_list', general_cuts)
        self.tree_mc.Draw('>>h_mc_ev_list', general_cuts)
        h_data_ev_list = gROOT.FindObject('h_data_ev_list')
        h_mc_ev_list = gROOT.FindObject('h_mc_ev_list')

        total_entries_data = h_data_ev_list.GetN()
        total_entries_mc = h_mc_ev_list.GetN()

        x_arr = np.arange(0, 5, 0.05)

        efficiency_data = np.zeros_like(x_arr)
        efficiency_mc = np.zeros_like(x_arr)

        for i, x in enumerate(x_arr):
            selection_cuts = 'cal_n_clusters == 1 && tr1_n_clusters == 1 && cal_cluster_y[0] > 158 && cal_cluster_y[0] < 168\
                    && abs(tr1_cluster_y - cal_cluster_y[0]) < 1.8 && abs(tr2_cluster_y - cal_cluster_y[0]) < {}'.format(x)

            self.tree_data.Draw('>>h_data_ev_list', selection_cuts)
            self.tree_mc.Draw('>>h_mc_ev_list', selection_cuts)
            h_data_ev_list = gROOT.FindObject('h_data_ev_list')
            h_mc_ev_list = gROOT.FindObject('h_mc_ev_list')

            efficiency_data[i] = h_data_ev_list.GetN() / total_entries_data * 100
            efficiency_mc[i] = h_mc_ev_list.GetN() / total_entries_mc * 100

        efficiency_data_graph = TGraph(len(x_arr), x_arr, efficiency_data)
        efficiency_data_graph.SetTitle('Data')
        efficiency_data_graph.Draw()

        efficiency_mc_graph = TGraph(len(x_arr), x_arr, efficiency_mc)
        efficiency_mc_graph.SetTitle('MC')
        efficiency_mc_graph.SetLineColor(2)
        efficiency_mc_graph.SetMarkerColor(2)
        efficiency_mc_graph.Draw('sameLP')

        canvas.BuildLegend()
        eff_data_graph_copy = efficiency_data_graph.Clone()
        eff_mc_graph_copy = efficiency_mc_graph.Clone()
        efficiency_data_graph.SetTitle('Tracker2 efficiency;|y_{tr2} - y_{shower}|;Efficiency, %')
        # Scaled subpad

        subpad = TPad("subpad", "", 0.35, 0.35, 0.9, 0.8)
        subpad.Draw()
        subpad.cd()
        subpad.SetGridx()
        subpad.SetGridy()

        eff_data_graph_copy.Draw()
        eff_data_graph_copy.SetMarkerStyle(7)
        eff_data_graph_copy.SetMaximum(100.5)
        eff_data_graph_copy.SetMinimum(98)
        eff_data_graph_copy.GetXaxis().SetLimits(1, 5)
        eff_data_graph_copy.SetTitle('')

        eff_mc_graph_copy.Draw('samePL')
        eff_mc_graph_copy.SetMarkerStyle(7)

        canvas.Write('tr2_efficiency')

    def backscattered_tracks(self):
        canvas = TCanvas('backscattered_tracks', 'title', 1024, 768)
        canvas_mc = TCanvas('backscattered_tracks_mc', 'title', 1024, 768)
        canvas.cd()
        gStyle.SetOptStat(1110)

        cuts = 'cal_n_clusters == 1 && cal_cluster_rho[0] > 158 && cal_cluster_rho[0] < 168'

        self.tree_data.Draw('(tr1_cluster_y[1]+23*4.5*(tr2_cluster_y[1]-tr1_cluster_y[1])/22.5):\
                             (tr1_cluster_x[1]+23*4.5*(tr2_cluster_x[1]-tr1_cluster_x[1])/22.5)\
                              >>h_tracks(100, -51, 51, 100, 0, 196)', cuts)

        h_tracks = gROOT.FindObject('h_tracks')
        h_tracks.SetTitle('Data, backscattered track hits in calorimeter;x, [mm];y, [mm]')
        h_tracks.Draw('colz')
        canvas.Write('h_tracks')

        canvas_mc.cd()
        self.tree_mc.Draw('(tr1_cluster_y[1]+23*4.5*(tr2_cluster_y[1]-tr1_cluster_y[1])/22.5):\
                             (tr1_cluster_x[1]+23*4.5*(tr2_cluster_x[1]-tr1_cluster_x[1])/22.5)\
                              >>h_tracks_mc(100, -51, 51, 100, 0, 196)', cuts)

        h_tracks_mc = gROOT.FindObject('h_tracks_mc')
        h_tracks_mc.SetTitle('MC, backscattered track hits in calorimeter;x, [mm];y, [mm]')

        h_tracks_mc.Draw('colz')
        canvas_mc.Write('h_tracks_mc')

    def tr1_distance_to_shower_scattered(self):
        canvas = TCanvas('tr1_distance_to_shower_scattered', 'title', 1024, 768)
        n_bins = 100
        first = -70
        last = 50

        self.tree_data.Draw('tr1_cluster_y[1]-cal_cluster_y[0]>>h_data({}, {}, {})'.format(n_bins, first, last))
        self.tree_mc.Draw('tr1_cluster_y[1]-cal_cluster_y[0]>>h_mc({}, {}, {})'.format(n_bins, first, last))
        h_data = gROOT.FindObject('h_data')
        h_mc = gROOT.FindObject('h_mc')

        h_mc.Draw('histo')
        h_mc.Scale(1. / h_mc.GetEntries())
        h_data.Draw('histosame')
        h_data.Scale(1. / h_data.GetEntries())

        h_mc.SetMaximum(max(h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin())))

        h_mc.SetFillColor(5)
        h_mc.SetTitle('MC;d_{cluster-shower}, [mm];#frac{N_{ev}}{N_{tot}}')

        h_data.SetLineWidth(3)
        h_data.SetTitle('Data')

        canvas.BuildLegend()
        h_mc.SetTitle('Distance:Tracker1 secondary cluster to the shower')
        gStyle.SetOptStat(1110)
        canvas.Write('tr1_distance_to_shower_scattered')

    def tr2_distance_to_shower_scattered(self):
        canvas = TCanvas('tr2_distance_to_shower_scattered', 'title', 1024, 768)
        n_bins = 100
        first = -70
        last = 50

        self.tree_data.Draw('tr2_cluster_y[1]-cal_cluster_y[0]>>h_data({}, {}, {})'.format(n_bins, first, last))
        self.tree_mc.Draw('tr2_cluster_y[1]-cal_cluster_y[0]>>h_mc({}, {}, {})'.format(n_bins, first, last))
        h_data = gROOT.FindObject('h_data')
        h_mc = gROOT.FindObject('h_mc')

        h_mc.Draw('histo')
        h_mc.Scale(1. / h_mc.GetEntries())
        h_data.Draw('histosame')
        h_data.Scale(1. / h_data.GetEntries())

        h_mc.SetMaximum(max(h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin())))

        h_mc.SetFillColor(5)
        h_mc.SetTitle('MC;d_{cluster-shower}, [mm];#frac{N_{ev}}{N_{tot}}')

        h_data.SetLineWidth(3)
        h_data.SetTitle('Data')

        canvas.BuildLegend()
        h_mc.SetTitle('Distance:Tracker2 secondary cluster to the shower')
        gStyle.SetOptStat(1110)
        canvas.Write('tr2_distance_to_shower_scattered')

    def all_clst_tr1_shower_distance(self):
        canvas = TCanvas('all_clst_tr1_shower_distance', 'title', 1024, 768)
        ratio_canvas = TCanvas('ratio_fits_tr1', 'title', 1024, 768)

        h_name = 'h_all_clst_tr1_shower_distance'

        cuts = 'cal_n_clusters == 1'
        distance = 'tr1_cluster_y - cal_cluster_y[0]'

        fit_peak = TF1("fit_peak", "gaus", -50, 30)
        fit_peak.SetMarkerStyle(1)
        fit_halo = TF1("fit_halo", "(0.5 * [0] * [1] / pi) / ((x[0] - [2]) * (x[0] - [2]) + .25 * [1] * [1])", -50, 30, 3)
        fit_halo.SetMarkerStyle(1)
        fit_total = TF1("fit_total", "gaus(0)+(0.5 * [3] * [4] / pi) / ((x[0] - [5]) * (x[0] - [5]) + .25 * [4] * [4])", -50, 30, 6)
        fit_total.SetMarkerStyle(1)

        func_ratio = TF1("func_ratio", "100 * (0.5 * [3] * [4] / pi) / ((x[0] - [5]) * (x[0] - [5]) + .25 * [4] * [4])/(gaus(0)+(0.5 * [3] * [4] / pi) / ((x[0] - [5]) * (x[0] - [5]) + .25 * [4] * [4]))", -50, 30, 6)
        func_ratio_mc = TF1("func_ratio_mc", "100 * (0.5 * [3] * [4] / pi) / ((x[0] - [5]) * (x[0] - [5]) + .25 * [4] * [4])/(gaus(0)+(0.5 * [3] * [4] / pi) / ((x[0] - [5]) * (x[0] - [5]) + .25 * [4] * [4]))", -50, 30, 6)

        fit_total.SetParameters(0.575, 0., 0.8, 0.06, 16., 0.6)

        fit_total.SetParLimits(0, 0.2, 0.9)
        fit_total.SetParLimits(1, -0.1, 0.1)
        fit_total.SetParLimits(2, 0.5, 1.5)
        fit_total.SetParLimits(3, 0.05, 0.2)
        fit_total.SetParLimits(4, 7, 20)
        fit_total.SetParLimits(5, -1.0, 1.0)
        fit_total.SetNpx(500)
        fit_peak.SetNpx(500)
        fit_halo.SetNpx(500)
        func_ratio.SetNpx(500)
        func_ratio_mc.SetNpx(500)

        self.tree_data.Draw(distance + '>>' + h_name + '(100, -70, 50)', cuts)
        self.tree_mc.Draw(distance + '>>' + h_name + '_mc' + '(100, -70, 50)', cuts)

        h_data = gROOT.FindObject(h_name)
        h_mc = gROOT.FindObject(h_name + '_mc')

        h_mc.Scale(1. / h_mc.GetEntries())
        h_data.Scale(1. / h_data.GetEntries())

        h_data.Fit(fit_total)

        par1 = array.array('d', 3 * [0.])
        par2 = array.array('d', 3 * [0.])
        pars = fit_total.GetParameters()

        par1[0], par1[1], par1[2] = pars[0], pars[1], pars[2]
        par2[0], par2[1], par2[2] = pars[3], pars[4], pars[5]

        fit_peak.SetParameters(par1)
        fit_halo.SetParameters(par2)
        func_ratio.SetParameters(pars)

        h_mc.Draw('histo')
        h_data.Draw('histosame')
        #fit_total.Draw('same')
        fit_peak.Draw('same')
        fit_halo.Draw('same')

        h_mc.SetMaximum(max(h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin())))
        h_mc.SetTitle('MC;d_{cluster-shower}, [mm];#frac{N_{ev}}{N_{tot}}')
        h_mc.SetFillColor(5)

        h_data.SetTitle('Data')
        h_data.SetLineWidth(3)

        fit_total.SetLineColor(2)
        fit_total.SetTitle('Total fit')
        fit_peak.SetLineColor(6)
        fit_peak.SetTitle('Peak fit')
        fit_halo.SetLineColor(7)
        fit_halo.SetTitle('Halo fit')

        canvas.SetLogy()
        canvas.BuildLegend()
        h_mc.SetTitle('Tracker1')

        canvas.Write(h_name)

        ratio_canvas.cd()
        h_mc.Fit(fit_total)
        pars_mc = fit_total.GetParameters()
        func_ratio_mc.SetParameters(pars_mc)

        func_ratio.Draw()
        func_ratio.SetTitle('Tracker1;d_{cluster-shower}, [mm];Probability, %')
        func_ratio_mc.Draw('same')
        func_ratio_mc.SetLineColor(2)
        ratio_canvas.Write('func_ratio_tr1')

    def all_clst_tr2_shower_distance(self):
        canvas = TCanvas('all_clst_tr2_shower_distance', 'title', 1024, 768)
        ratio_canvas = TCanvas('ratio_fits_tr2', 'title', 1024, 768)

        h_name = 'h_all_clst_tr2_shower_distance'

        cuts = 'cal_n_clusters == 1'
        distance = 'tr2_cluster_y - cal_cluster_y[0]'

        fit_peak = TF1("fit_peak", "gaus", -50, 30)
        fit_peak.SetMarkerStyle(1)
        fit_halo = TF1("fit_halo", "(0.5 * [0] * [1] / pi) / ((x[0] - [2]) * (x[0] - [2]) + .25 * [1] * [1])", -50, 30, 3)
        fit_halo.SetMarkerStyle(1)
        fit_total = TF1("fit_total", "gaus(0)+(0.5 * [3] * [4] / pi) / ((x[0] - [5]) * (x[0] - [5]) + .25 * [4] * [4])", -50, 30, 6)
        fit_total.SetMarkerStyle(1)

        func_ratio = TF1("func_ratio", "100 * (0.5 * [3] * [4] / pi) / ((x[0] - [5]) * (x[0] - [5]) + .25 * [4] * [4])/(gaus(0)+(0.5 * [3] * [4] / pi) / ((x[0] - [5]) * (x[0] - [5]) + .25 * [4] * [4]))", -50, 30, 6)
        func_ratio_mc = TF1("func_ratio_mc", "100 * (0.5 * [3] * [4] / pi) / ((x[0] - [5]) * (x[0] - [5]) + .25 * [4] * [4])/(gaus(0)+(0.5 * [3] * [4] / pi) / ((x[0] - [5]) * (x[0] - [5]) + .25 * [4] * [4]))", -50, 30, 6)

        fit_total.SetParameters(0.575, 0., 0.8, 0.06, 16., 0.6)

        fit_total.SetParLimits(0, 0.2, 0.9)
        fit_total.SetParLimits(1, -0.1, 0.1)
        fit_total.SetParLimits(2, 0.5, 1.5)
        fit_total.SetParLimits(3, 0.05, 0.2)
        fit_total.SetParLimits(4, 7, 20)
        fit_total.SetParLimits(5, -1.0, 1.0)
        fit_total.SetNpx(500)
        fit_peak.SetNpx(500)
        fit_halo.SetNpx(500)
        func_ratio.SetNpx(500)
        func_ratio_mc.SetNpx(500)

        self.tree_data.Draw(distance + '>>' + h_name + '(100, -70, 50)', cuts)
        self.tree_mc.Draw(distance + '>>' + h_name + '_mc' + '(100, -70, 50)', cuts)

        h_data = gROOT.FindObject(h_name)
        h_mc = gROOT.FindObject(h_name + '_mc')

        h_mc.Scale(1. / h_mc.GetEntries())
        h_data.Scale(1. / h_data.GetEntries())

        h_data.Fit(fit_total)

        par1 = array.array('d', 3 * [0.])
        par2 = array.array('d', 3 * [0.])
        pars = fit_total.GetParameters()

        par1[0], par1[1], par1[2] = pars[0], pars[1], pars[2]
        par2[0], par2[1], par2[2] = pars[3], pars[4], pars[5]

        fit_peak.SetParameters(par1)
        fit_halo.SetParameters(par2)
        func_ratio.SetParameters(pars)

        h_mc.Draw('histo')
        h_data.Draw('histosame')
        #fit_total.Draw('same')
        fit_peak.Draw('same')
        fit_halo.Draw('same')

        h_mc.SetMaximum(max(h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin())))
        h_mc.SetTitle('MC;d_{cluster-shower}, [mm];#frac{N_{ev}}{N_{tot}}')
        h_mc.SetFillColor(5)

        h_data.SetTitle('Data')
        h_data.SetLineWidth(3)

        fit_total.SetLineColor(2)
        fit_total.SetTitle('Total fit')
        fit_peak.SetLineColor(6)
        fit_peak.SetTitle('Peak fit')
        fit_halo.SetLineColor(7)
        fit_halo.SetTitle('Halo fit')

        canvas.SetLogy()
        canvas.BuildLegend()
        h_mc.SetTitle('Tracker2')

        canvas.Write(h_name)
        ratio_canvas.cd()

        h_mc.Fit(fit_total)
        pars_mc = fit_total.GetParameters()
        func_ratio_mc.SetParameters(pars_mc)

        func_ratio.Draw()
        func_ratio.SetTitle('Tracker2;d_{cluster-shower}, [mm];Probability, %')
        func_ratio_mc.Draw('same')
        func_ratio_mc.SetLineColor(2)
        ratio_canvas.Write('func_ratio_tr2')

    def simple_event(self):
        canvas = TCanvas('simple_event', 'title', 1024, 768)
        canvas.cd()
        gStyle.SetOptStat(0)

        cuts_event = 'cal_tower_energy*(Entry$==532)'
        cuts_cluster = '(nomerge.cal_tower_cluster+1)*(Entry$==532)'
        cuts_cluster_merged = '(cal_tower_cluster+1)*(Entry$==532)'

        self.tree_data.Draw('cal_tower_pad:cal_tower_sector>>h_event(5,0,5,20,35,55)', cuts_event)
        h_event = gROOT.FindObject('h_event')
        h_event.Draw('colztext')
        h_event.SetTitle('Signals, event #532;sector number;pad number')
        canvas.Write('h_event')

        self.tree_data.Draw('nomerge.cal_tower_pad:nomerge.cal_tower_sector>>h_cluster(5,0,5,20,35,55)', cuts_cluster)
        h_cluster = gROOT.FindObject('h_cluster')
        h_cluster.Draw('colztext')
        h_cluster.SetTitle('Clusters, event #532;sector number;pad number')
        canvas.Write('h_cluster')

        self.tree_data.Draw('cal_tower_pad:cal_tower_sector>>h_cluster_merged(5,0,5,20,35,55)', cuts_cluster_merged)
        h_cluster_merged = gROOT.FindObject('h_cluster_merged')
        h_cluster_merged.Draw('colztext')
        h_cluster_merged.SetTitle('Clusters merged, event #532;sector number;pad number')
        canvas.Write('h_cluster_merged')

    def number_of_bs_events(self):
        good_cut = 'cal_n_clusters == 1'
        cuts_11 = 'tr1_n_clusters == 1 && tr2_n_clusters == 1 && cal_n_clusters == 1'
        cuts_12 = 'tr1_n_clusters == 1 && tr2_n_clusters == 2 && cal_n_clusters == 1'
        cuts_21 = 'tr1_n_clusters == 2 && tr2_n_clusters == 1 && cal_n_clusters == 1'
        cuts_22 = 'tr1_n_clusters == 2 && tr2_n_clusters == 2 && cal_n_clusters == 1'
        cuts_10 = 'tr1_n_clusters == 1 && tr2_n_clusters == 0 && cal_n_clusters == 1'
        cuts_01 = 'tr1_n_clusters == 0 && tr2_n_clusters == 1 && cal_n_clusters == 1'
        cuts_00 = 'tr1_n_clusters == 0 && tr2_n_clusters == 0 && cal_n_clusters == 1'

        self.tree_data.Draw('>>good_data', good_cut)
        self.tree_data.Draw('>>h11_data', cuts_11)
        self.tree_data.Draw('>>h12_data', cuts_12)
        self.tree_data.Draw('>>h21_data', cuts_21)
        self.tree_data.Draw('>>h22_data', cuts_22)
        self.tree_data.Draw('>>h10_data', cuts_10)
        self.tree_data.Draw('>>h01_data', cuts_01)
        self.tree_data.Draw('>>h00_data', cuts_00)

        n_data_events = gROOT.FindObject('good_data').GetN()
        print('Total DATA Events:', n_data_events)

        print('11 Events:', 100 * gROOT.FindObject('h11_data').GetN() / n_data_events)
        print('12 Events:', 100 * gROOT.FindObject('h12_data').GetN() / n_data_events)
        print('21 Events:', 100 * gROOT.FindObject('h21_data').GetN() / n_data_events)
        print('22 Events:', 100 * gROOT.FindObject('h22_data').GetN() / n_data_events)
        print('10 Events:', 100 * gROOT.FindObject('h10_data').GetN() / n_data_events)
        print('01 Events:', 100 * gROOT.FindObject('h01_data').GetN() / n_data_events)
        print('00 Events:', 100 * gROOT.FindObject('h00_data').GetN() / n_data_events)

        self.tree_mc.Draw('>>good_mc', good_cut)
        self.tree_mc.Draw('>>h11_mc', cuts_11)
        self.tree_mc.Draw('>>h12_mc', cuts_12)
        self.tree_mc.Draw('>>h21_mc', cuts_21)
        self.tree_mc.Draw('>>h22_mc', cuts_22)
        self.tree_mc.Draw('>>h10_mc', cuts_10)
        self.tree_mc.Draw('>>h01_mc', cuts_01)
        self.tree_mc.Draw('>>h00_mc', cuts_00)
        n_mc_events = gROOT.FindObject('good_mc').GetN()

        print('Total mc Events:', n_mc_events)
        print('11 Events:', 100 * gROOT.FindObject('h11_mc').GetN() / n_mc_events)
        print('12 Events:', 100 * gROOT.FindObject('h12_mc').GetN() / n_mc_events)
        print('21 Events:', 100 * gROOT.FindObject('h21_mc').GetN() / n_mc_events)
        print('22 Events:', 100 * gROOT.FindObject('h22_mc').GetN() / n_mc_events)
        print('10 Events:', 100 * gROOT.FindObject('h10_mc').GetN() / n_mc_events)
        print('01 Events:', 100 * gROOT.FindObject('h01_mc').GetN() / n_mc_events)
        print('00 Events:', 100 * gROOT.FindObject('h00_mc').GetN() / n_mc_events)

    def w0_resolution(self):
        canvas = TCanvas('w0_scan', 'title', 1024, 768)
        canvas.cd()
        n = 9
        x = np.array([3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8])
        y = np.array([1.72355, 1.71761, 1.71620, 1.71436, 1.71386, 1.71230, 1.71305, 1.71492, 1.72129])
        w0_graph = TGraph(n, x, y)
        w0_graph.Draw()
        canvas.Write('w0_graph')


    def shower_distance(self, tr_number, cluster_number, tr1_n_clusters, tr2_n_clusters):
        canvas = TCanvas('shower_distance', 'title', 1024, 768)

        h_name = 'tr{}_clst{}_shower_distance_{}{}'.format(tr_number, cluster_number, tr1_n_clusters, tr2_n_clusters)

        distance = 'tr{}_cluster_rho[{}] - cal_cluster_rho[0]'.format(tr_number, cluster_number - 1)

        cuts = 'tr1_n_clusters == {} && tr2_n_clusters == {}'.format(tr1_n_clusters, tr2_n_clusters)

        self.tree_data.Draw(distance + '>>' + h_name + '(100, -60, 40)', cuts)
        self.tree_mc.Draw(distance + '>>' + h_name + '_mc' + '(100, -60, 40)', cuts)

        h_data = gROOT.FindObject(h_name)
        h_mc = gROOT.FindObject(h_name + '_mc')

        h_mc.Draw('histo')
        h_mc.Scale(1. / self.n_events_mc)
        h_data.Draw('histosame')
        h_data.Scale(1. / self.n_events_data)

        h_mc.SetMaximum(max(h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin())))
        h_mc.SetTitle('MC')
        h_mc.GetXaxis().SetTitle('Distance to the shower, mm')
        h_mc.GetYaxis().SetTitle('N_{events} normalized')
        h_mc.SetFillColor(5)

        h_data.SetTitle('Data')
        h_data.SetLineWidth(3)

        canvas.SetLogy()

        canvas.BuildLegend()
        h_mc.SetTitle(h_name)

        canvas.Write(h_name)

    def residuals_111(self):
        z_tr1 = 0.
        z_tr2 = 5 * 4.5
        z_cal = 25 * 4.5
        z = np.array([z_tr1, z_tr2, z_cal])
        z_err = np.array([4.5 / 2, 4.5 / 2, 13.5])
        rho_err = np.array([0.9, 0.9, 0.9])

        canvas = TCanvas('residuals_111', 'title', 1024, 1024)

        h_tr1 = TH1F('h_tr1_res111', 'Tracker1 residuals', 500, -20, 20)
        h_tr2 = TH1F('h_tr2_res111', 'Tracker2 residuals', 500, -20, 20)
        h_cal = TH1F('h_cal_res111', 'Calorimeter residuals', 500, -20, 20)

        h_tr1_mc = TH1F('h_tr1_res111_mc', 'Tracker1 residuals MC', 500, -20, 20)
        h_tr2_mc = TH1F('h_tr2_res111_mc', 'Tracker2 residuals MC', 500, -20, 20)
        h_cal_mc = TH1F('h_cal_res111_mc', 'Calorimeter residuals MC', 500, -20, 20)

        for event in self.tree_data:
            if not (event.tr1_n_clusters == 1
                    and event.tr2_n_clusters == 1
                    and event.cal_n_clusters == 1
                    and 153.1 < event.cal_cluster_rho[0] < 172.3):
                continue

            rho = np.array([event.tr1_cluster_rho[0], event.tr2_cluster_rho[0], event.cal_cluster_rho[0]])
            track = TGraphErrors(3, z, rho, z_err, rho_err)
            track.Fit('pol0', "Q")
            fit_func = track.GetFunction('pol0')

            tr1_residual = event.tr1_cluster_rho[0] - fit_func.Eval(z_tr1)
            tr2_residual = event.tr2_cluster_rho[0] - fit_func.Eval(z_tr2)
            cal_residual = event.cal_cluster_rho[0] - fit_func.Eval(z_cal)
            h_tr1.Fill(tr1_residual)
            h_tr2.Fill(tr2_residual)
            h_cal.Fill(cal_residual)

        for event in self.tree_mc:
            if not (event.tr1_n_clusters == 1
                    and event.tr2_n_clusters == 1
                    and event.cal_n_clusters == 1
                    and 153.1 < event.cal_cluster_rho[0] < 172.3):
                continue

            rho = np.array([event.tr1_cluster_rho[0], event.tr2_cluster_rho[0], event.cal_cluster_rho[0]])
            track = TGraphErrors(3, z, rho, z_err, rho_err)
            track.Fit('pol0', "Q")
            fit_func = track.GetFunction('pol0')

            tr1_residual = event.tr1_cluster_rho[0] - fit_func.Eval(z_tr1)
            tr2_residual = event.tr2_cluster_rho[0] - fit_func.Eval(z_tr2)
            cal_residual = event.cal_cluster_rho[0] - fit_func.Eval(z_cal)
            h_tr1_mc.Fill(tr1_residual)
            h_tr2_mc.Fill(tr2_residual)
            h_cal_mc.Fill(cal_residual)

        h_tr1_mc.Draw('histo')
        h_tr1_mc.Scale(1. / self.n_events_mc)

        h_tr1_mc.SetFillColor(5)

        h_tr1.Draw('histosame')
        h_tr1.Scale(1. / self.n_events_data)
        h_tr1.SetLineWidth(3)

        h_tr1_mc.SetMaximum(max(h_tr1_mc.GetBinContent(h_tr1_mc.GetMaximumBin()), h_tr1.GetBinContent(h_tr1.GetMaximumBin())))
        canvas.BuildLegend()
        canvas.Write('tr1_res111')

        h_tr2_mc.Draw('histo')
        h_tr2_mc.Scale(1. / self.n_events_mc)

        h_tr2_mc.SetFillColor(5)

        h_tr2.Draw('histosame')
        h_tr2.Scale(1. / self.n_events_data)
        h_tr2.SetLineWidth(3)

        h_tr2_mc.SetMaximum(max(h_tr2_mc.GetBinContent(h_tr2_mc.GetMaximumBin()), h_tr2.GetBinContent(h_tr2.GetMaximumBin())))
        canvas.BuildLegend()
        canvas.Write('tr2_res111')

        h_cal_mc.Draw('histo')
        h_cal_mc.Scale(1. / self.n_events_mc)

        h_cal_mc.SetFillColor(5)

        h_cal.Draw('histosame')
        h_cal.Scale(1. / self.n_events_data)
        h_cal.SetLineWidth(3)

        h_cal_mc.SetMaximum(max(h_cal_mc.GetBinContent(h_cal_mc.GetMaximumBin()), h_cal.GetBinContent(h_cal.GetMaximumBin())))
        canvas.BuildLegend()
        canvas.Write('cal_res111')

    def distance_between_clusters(self, det, clst1, clst2):
        canvas = TCanvas('distance_between_clusters', 'distance between clusters', 1024, 1024)
        n_bins = 128
        first = -115.2 / 2
        last = 115.2 / 2

        distance = det + '_cluster_rho[{}] -'.format(clst1) + det + '_cluster_rho[{}]'.format(clst2)

        self.tree_data.Draw(distance + '>>h_data({}, {}, {})'.format(n_bins, first, last))  # 'cal_n_clusters == 1'
        self.tree_mc.Draw(distance + '>>h_mc({}, {}, {})'.format(n_bins, first, last))  # 'cal_n_clusters == 1'
        h_mc = gROOT.FindObject('h_mc')
        h_data = gROOT.FindObject('h_data')

        h_mc.Draw('histo')
        h_mc.Scale(1. / self.n_events_mc)
        h_data.Draw('histosame')
        h_data.Scale(1. / self.n_events_data)
        h_mc.SetMaximum(max(h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin())))

        h_mc.SetFillColor(5)
        h_mc.SetTitle('MC')

        h_data.SetLineWidth(3)
        h_data.SetTitle('Data')

        canvas.BuildLegend()
        canvas.Write('clst{}_clst{}_distance in '.format(clst1, clst2) + det)

    def resolution(self, det):
        n_bins = 120
        first = -30
        last = 30

        cuts = 'cal_n_clusters == 1 && cal_cluster_rho[0] > 153.1 && cal_cluster_rho[0] < 172.3 && cal_cluster_energy[0] > 100'
        self.tree_mc.Draw(det + '_cluster_rho[0] - ' + det + '_true_hit_rho>>h_mc({}, {}, {})'.format(n_bins, first, last), cuts)
        h_mc = gROOT.FindObject('h_mc')

        h_mc.Draw('histo')
        h_mc.Scale(1. / self.n_events_mc)

        h_mc.SetFillColor(5)
        h_mc.SetTitle('MC')
        h_mc.Write('resolution_of_' + det)

    def n_clusters(self, det):
        canvas = TCanvas('n_clusters', 'title', 1024, 768)
        n_bins = 15
        first = 0
        last = 15

        cuts = 'cal_n_clusters == 1 && cal_cluster_rho[0] > 153.1 && cal_cluster_rho[0] < 172.3 && cal_cluster_energy[0] > 100'

        self.tree_data.Draw(det + '_n_clusters>>h_data({}, {}, {})'.format(n_bins, first, last), cuts)
        self.tree_mc.Draw(det + '_n_clusters>>h_mc({}, {}, {})'.format(n_bins, first, last), cuts)

        h_data = gROOT.FindObject('h_data')
        h_mc = gROOT.FindObject('h_mc')

        h_mc.Draw('histo')
        h_mc.Scale(1. / self.n_events_mc)
        h_data.Draw('histosame')
        h_data.Scale(1. / self.n_events_data)

        h_mc.SetMaximum(max(h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin())))
        h_mc.SetTitle('MC')
        h_mc.GetXaxis().SetTitle('N clusters')
        h_mc.GetYaxis().SetTitle('N_{events} normalized')
        h_mc.SetFillColor(5)

        h_data.SetTitle('Data')
        h_data.SetLineWidth(3)

        canvas.BuildLegend()

        canvas.Write(det + '_n_clusters')

    def clusters_energy(self):
        for idx, tree in enumerate([self.tree_data, self.tree_mc]):
            canvas = TCanvas('clusters_energy' + str(idx), 'title', 1024, 768)
            n_bins = 500
            first = 0
            last = 500

            tree.Draw('cal_cluster_energy[0]>>h_clst1_energy({}, {}, {})'.format(n_bins, first, last))
            tree.Draw('cal_cluster_energy[1]>>h_clst2_energy({}, {}, {})'.format(n_bins, first, last))
            tree.Draw('cal_cluster_energy[2]>>h_clst3_energy({}, {}, {})'.format(n_bins, first, last))

            h_clst1_energy = gROOT.FindObject('h_clst1_energy')
            h_clst2_energy = gROOT.FindObject('h_clst2_energy')
            h_clst3_energy = gROOT.FindObject('h_clst3_energy')

            h_clst1_energy.Scale(1. / self.n_events_data)
            h_clst2_energy.Scale(1. / self.n_events_data)
            h_clst3_energy.Scale(1. / self.n_events_data)

            h_clst1_energy.Draw('histo')
            h_clst2_energy.Draw('histosame')
            h_clst3_energy.Draw('histosame')

            h_clst1_energy.SetMaximum(max([h_clst1_energy.GetBinContent(h_clst1_energy.GetMaximumBin()), h_clst2_energy.GetBinContent(h_clst2_energy.GetMaximumBin()), h_clst3_energy.GetBinContent(h_clst3_energy.GetMaximumBin())]))

            h_clst1_energy.SetTitle('Cluster1')
            h_clst2_energy.SetTitle('Cluster2')
            h_clst3_energy.SetTitle('Cluster3')

            h_clst1_energy.GetXaxis().SetTitle('Energy, [MIP]')
            h_clst1_energy.GetYaxis().SetTitle('N_{events} normalized')
            h_clst1_energy.SetLineColor(1)
            h_clst2_energy.SetLineColor(2)
            h_clst3_energy.SetLineColor(3)
            h_clst1_energy.SetLineWidth(3)
            h_clst2_energy.SetLineWidth(3)
            h_clst3_energy.SetLineWidth(3)

            canvas.BuildLegend()
            if idx == 0:
                canvas.Write('clusters_energies_data')
            else:
                canvas.Write('clusters_energies_MC')

    def merge_pos_change(self):
        canvas = TCanvas('merge_pos_change', 'title', 1024, 768)

        self.tree_data.Draw('nomerge.cal_cluster_rho[0] - cal_cluster_rho[0]>>h_data(300,0,6)')
        self.tree_mc.Draw('nomerge.cal_cluster_rho[0] - cal_cluster_rho[0]>>h_mc(300,0,6)')

        h_data = gROOT.FindObject('h_data')
        h_mc = gROOT.FindObject('h_mc')

        h_mc.Draw('histo')
        h_mc.Scale(1. / self.n_events_mc)
        h_data.Draw('histosame')
        h_data.Scale(1. / self.n_events_data)

        h_mc.SetMaximum(max(h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin())))
        h_mc.SetTitle('MC')
        h_mc.GetXaxis().SetTitle('position shift, [mm]')
        h_mc.GetYaxis().SetTitle('N_{events} normalized')
        h_mc.SetFillColor(5)

        h_data.SetTitle('Data')
        h_data.SetLineWidth(3)

        canvas.BuildLegend()
        canvas.Write('merge_pos_change')


def shower_distances(data):
    '''
    Calculate efficiency for Tr2 if Tr1 has signal
    in range of +-"sigma" = 4.4 mm within shower in cal
    '''

    # tr_number, cluster_number, tr1_n_clusters, tr2_n_clusters
    # 111
    data.shower_distance(1, 1, 1, 1)
    data.shower_distance(2, 1, 1, 1)
    # 121
    data.shower_distance(1, 1, 1, 2)
    data.shower_distance(2, 1, 1, 2)
    data.shower_distance(2, 2, 1, 2)
    # 211
    data.shower_distance(1, 1, 2, 1)
    data.shower_distance(1, 2, 2, 1)
    data.shower_distance(2, 1, 2, 1)
    # 221
    data.shower_distance(1, 1, 2, 2)
    data.shower_distance(1, 2, 2, 2)
    data.shower_distance(2, 1, 2, 2)
    data.shower_distance(2, 2, 2, 2)


def main():
    # data.simple_event()
    # data.number_of_bs_events()
    # data.w0_resolution()
    data = Data()
    data.calorimeter_plot_x()
    data.calorimeter_plot_y()
    data.calorimeter_plot_layer()
    data.calorimeter_plot_n_clusters()
    data.calorimeter_plot_n_pads()
    data.calorimeter_plot_energy()
    data.tr1_distance_to_shower()
    data.tr2_distance_to_shower()
    data.tr1_efficiency()
    data.tr2_efficiency()
    data.backscattered_tracks()
    data.tr1_distance_to_shower_scattered()
    data.tr2_distance_to_shower_scattered()
    data.all_clst_tr1_shower_distance()
    data.all_clst_tr2_shower_distance()


    # data.calorimeter_plot_n_pads_vs_energy()

    input('Yaay I am finished :3')


gROOT.SetBatch(1)
gROOT.SetStyle('ATLAS')
gStyle.SetOptStat(1110)

main()
