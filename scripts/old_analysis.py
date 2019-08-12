from ROOT import TFile, gROOT, TGraphErrors, TH1F, TGraph, TCanvas, TPad, TF1, gStyle
import numpy as np
import array
from store_file import langaufun


class Detector:
    file_data = TFile.Open('./extracted_trees/with_1.4_energy_cut/extracted_5_gev_energy_scan_1.root', 'read')
    tree_data = file_data.data

    # file_mc = TFile.Open('./result_trees/extracted_mc.root', 'read')
    # tree_mc = file_mc.mc

    output_file = TFile('histos.root', 'recreate')

    def __init__(self):
        self.output_file.cd()

    @classmethod
    def check_alignment(cls):
        x = np.array([0., 5 * 4.5, 25 * 4.5])

        cuts = "tr1_n_clusters == 1 && tr2_n_clusters == 1 && cal_n_clusters == 1 && 153.1 < cal_cluster_y[0] && cal_cluster_y[0] < 172.3"

        cls.tree_data.Draw('tr1_cluster_y[0]>>h_tr1(64, 80, 195.2)', cuts)
        h_tr1 = gROOT.FindObject('h_tr1')
        y_tr1 = h_tr1.GetMean()

        cls.tree_data.Draw('tr2_cluster_y[0]>>h_tr2(64, 80, 195.2)', cuts)
        h_tr2 = gROOT.FindObject('h_tr2')
        y_tr2 = h_tr2.GetMean()

        cls.tree_data.Draw('cal_cluster_y[0]>>h_cal(64, 80, 195.2)', cuts)
        h_cal = gROOT.FindObject('h_cal')
        y_cal = h_cal.GetMean()

        track = TGraphErrors(3, x, np.array([y_tr1, y_tr2, y_cal]))
        track.Fit('pol0', "Q")
        fit_func = track.GetFunction('pol0')
        print("Tr1 needs alignment for:", y_tr1 - fit_func.Eval(0))
        print("Tr2 needs alignment for:", y_tr2 - fit_func.Eval(0))
        print("Cal needs alignment for:", y_cal - fit_func.Eval(0))

    @classmethod
    def check_layer2_energy(cls):
        canvas = TCanvas('energy_layer2', 'energy_layer2', 1024, 768)
        n_bins = 100
        first = 0
        last = 200

        cls.tree_data.Draw('Sum$(cal_hit_energy*(cal_hit_layer == 3))>>h_data({}, {}, {})'.format(n_bins, first, last))
        cls.tree_mc.Draw('Sum$(cal_hit_energy*(cal_hit_layer == 3))>>h_mc({}, {}, {})'.format(n_bins, first, last))
        h_data = gROOT.FindObject('h_data')
        h_mc = gROOT.FindObject('h_mc')

        h_mc.Draw('histo')
        h_mc.Scale(1. / h_mc.GetEntries())
        h_data.Draw('histosame')
        h_data.Scale(1. / h_data.GetEntries())

        h_mc.SetMaximum(max(h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin())) + 0.01)

        h_mc.SetFillColor(5)
        h_mc.SetTitle('MC;energy, [MIP];#frac{N_{ev}}{N_{tot}}')

        h_data.SetLineWidth(3)
        h_data.SetTitle('Data')

        canvas.BuildLegend()
        h_mc.SetTitle('energy_layer2')
        gStyle.SetOptStat(1110)
        canvas.Write('energy_layer2')

    @classmethod
    def backscattered_tracks(cls):
        canvas = TCanvas('backscattered_tracks', 'backscattered_tracks', 1024, 768)
        canvas_mc = TCanvas('backscattered_tracks_mc', 'backscattered_tracks_mc', 1024, 768)
        canvas.cd()
        gStyle.SetOptStat(1110)

        cuts = 'cal_n_clusters == 1 && cal_cluster_rho[0] > 158 && cal_cluster_rho[0] < 168'

        cls.tree_data.Draw('(tr1_cluster_y[1]+23*4.5*(tr2_cluster_y[1]-tr1_cluster_y[1])/22.5):\
                             (tr1_cluster_x[1]+23*4.5*(tr2_cluster_x[1]-tr1_cluster_x[1])/22.5)\
                              >>h_tracks(100, -51, 51, 100, 0, 196)', cuts)

        h_tracks = gROOT.FindObject('h_tracks')
        h_tracks.SetTitle('Data, backscattered track hits in calorimeter;x, [mm];y, [mm]')
        h_tracks.Draw('colz')
        canvas.Write('backscattered_tracks')

        canvas_mc.cd()
        cls.tree_mc.Draw('(tr1_cluster_y[1]+23*4.5*(tr2_cluster_y[1]-tr1_cluster_y[1])/22.5):\
                             (tr1_cluster_x[1]+23*4.5*(tr2_cluster_x[1]-tr1_cluster_x[1])/22.5)\
                              >>h_tracks_mc(100, -51, 51, 100, 0, 196)', cuts)

        h_tracks_mc = gROOT.FindObject('h_tracks_mc')
        h_tracks_mc.SetTitle('MC, backscattered track hits in calorimeter;x, [mm];y, [mm]')

        h_tracks_mc.Draw('colz')
        canvas_mc.Write('backscattered_tracks_mc')

    @classmethod
    def simple_event(cls):
        canvas = TCanvas('simple_event', 'title', 1024, 768)
        canvas.cd()
        gStyle.SetOptStat(0)

        cuts_event = 'cal_tower_energy*(Entry$==532)'
        cuts_cluster = '(nomerge.cal_tower_cluster+1)*(Entry$==532)'
        cuts_cluster_merged = '(cal_tower_cluster+1)*(Entry$==532)'

        cls.tree_data.Draw('cal_tower_pad:cal_tower_sector>>h_event(5,0,5,20,35,55)', cuts_event)
        h_event = gROOT.FindObject('h_event')
        h_event.Draw('colztext')
        h_event.SetTitle('Signals, event #532;sector number;pad number')
        canvas.Write('h_event')

        cls.tree_data.Draw('nomerge.cal_tower_pad:nomerge.cal_tower_sector>>h_cluster(5,0,5,20,35,55)', cuts_cluster)
        h_cluster = gROOT.FindObject('h_cluster')
        h_cluster.Draw('colztext')
        h_cluster.SetTitle('Clusters, event #532;sector number;pad number')
        canvas.Write('h_cluster')

        cls.tree_data.Draw('cal_tower_pad:cal_tower_sector>>h_cluster_merged(5,0,5,20,35,55)', cuts_cluster_merged)
        h_cluster_merged = gROOT.FindObject('h_cluster_merged')
        h_cluster_merged.Draw('colztext')
        h_cluster_merged.SetTitle('Clusters merged, event #532;sector number;pad number')
        canvas.Write('h_cluster_merged')

    @classmethod
    def number_of_bs_events(cls):
        good_cut = 'cal_n_clusters == 1'
        cuts_11 = 'tr1_n_clusters == 1 && tr2_n_clusters == 1 && cal_n_clusters == 1'
        cuts_12 = 'tr1_n_clusters == 1 && tr2_n_clusters == 2 && cal_n_clusters == 1'
        cuts_21 = 'tr1_n_clusters == 2 && tr2_n_clusters == 1 && cal_n_clusters == 1'
        cuts_22 = 'tr1_n_clusters == 2 && tr2_n_clusters == 2 && cal_n_clusters == 1'
        cuts_10 = 'tr1_n_clusters == 1 && tr2_n_clusters == 0 && cal_n_clusters == 1'
        cuts_01 = 'tr1_n_clusters == 0 && tr2_n_clusters == 1 && cal_n_clusters == 1'
        cuts_00 = 'tr1_n_clusters == 0 && tr2_n_clusters == 0 && cal_n_clusters == 1'

        cls.tree_data.Draw('>>good_data', good_cut)
        cls.tree_data.Draw('>>h11_data', cuts_11)
        cls.tree_data.Draw('>>h12_data', cuts_12)
        cls.tree_data.Draw('>>h21_data', cuts_21)
        cls.tree_data.Draw('>>h22_data', cuts_22)
        cls.tree_data.Draw('>>h10_data', cuts_10)
        cls.tree_data.Draw('>>h01_data', cuts_01)
        cls.tree_data.Draw('>>h00_data', cuts_00)

        n_data_events = gROOT.FindObject('good_data').GetN()
        print('Total DATA Events:', n_data_events)

        print('11 Events:', 100 * gROOT.FindObject('h11_data').GetN() / n_data_events)
        print('12 Events:', 100 * gROOT.FindObject('h12_data').GetN() / n_data_events)
        print('21 Events:', 100 * gROOT.FindObject('h21_data').GetN() / n_data_events)
        print('22 Events:', 100 * gROOT.FindObject('h22_data').GetN() / n_data_events)
        print('10 Events:', 100 * gROOT.FindObject('h10_data').GetN() / n_data_events)
        print('01 Events:', 100 * gROOT.FindObject('h01_data').GetN() / n_data_events)
        print('00 Events:', 100 * gROOT.FindObject('h00_data').GetN() / n_data_events)

        cls.tree_mc.Draw('>>good_mc', good_cut)
        cls.tree_mc.Draw('>>h11_mc', cuts_11)
        cls.tree_mc.Draw('>>h12_mc', cuts_12)
        cls.tree_mc.Draw('>>h21_mc', cuts_21)
        cls.tree_mc.Draw('>>h22_mc', cuts_22)
        cls.tree_mc.Draw('>>h10_mc', cuts_10)
        cls.tree_mc.Draw('>>h01_mc', cuts_01)
        cls.tree_mc.Draw('>>h00_mc', cuts_00)
        n_mc_events = gROOT.FindObject('good_mc').GetN()

        print('Total mc Events:', n_mc_events)
        print('11 Events:', 100 * gROOT.FindObject('h11_mc').GetN() / n_mc_events)
        print('12 Events:', 100 * gROOT.FindObject('h12_mc').GetN() / n_mc_events)
        print('21 Events:', 100 * gROOT.FindObject('h21_mc').GetN() / n_mc_events)
        print('22 Events:', 100 * gROOT.FindObject('h22_mc').GetN() / n_mc_events)
        print('10 Events:', 100 * gROOT.FindObject('h10_mc').GetN() / n_mc_events)
        print('01 Events:', 100 * gROOT.FindObject('h01_mc').GetN() / n_mc_events)
        print('00 Events:', 100 * gROOT.FindObject('h00_mc').GetN() / n_mc_events)

    @staticmethod
    def w0_resolution():
        canvas = TCanvas('w0_scan', 'title', 1024, 768)
        canvas.cd()
        n = 9
        x = np.array([3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8])
        y = np.array([1.72355, 1.71761, 1.71620, 1.71436, 1.71386, 1.71230, 1.71305, 1.71492, 1.72129])
        w0_graph = TGraph(n, x, y)
        w0_graph.Draw()
        canvas.Write('w0_graph')

    @classmethod
    def shower_distance(cls, tr_number, cluster_number, tr1_n_clusters, tr2_n_clusters):
        canvas = TCanvas('shower_distance', 'title', 1024, 768)

        h_name = 'tr{}_clst{}_shower_distance_{}{}'.format(tr_number, cluster_number, tr1_n_clusters, tr2_n_clusters)

        distance = 'tr{}_cluster_rho[{}] - cal_cluster_rho[0]'.format(tr_number, cluster_number - 1)

        cuts = 'tr1_n_clusters == {} && tr2_n_clusters == {}'.format(tr1_n_clusters, tr2_n_clusters)

        cls.tree_data.Draw(distance + '>>' + h_name + '(100, -60, 40)', cuts)
        cls.tree_mc.Draw(distance + '>>' + h_name + '_mc' + '(100, -60, 40)', cuts)

        h_data = gROOT.FindObject(h_name)
        h_mc = gROOT.FindObject(h_name + '_mc')

        h_mc.Draw('histo')
        h_mc.Scale(1. / h_mc.GetEntries())
        h_data.Draw('histosame')
        h_data.Scale(1. / h_mc.GetEntries())

        h_mc.SetMaximum(max(h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin())) + 0.01)
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

    @classmethod
    def e_ratio_tr2_tr1(cls):
        canvas = TCanvas('e_ratio_tr2_tr1', 'e_ratio_tr2_tr1', 1024, 768)
        n_bins = 15
        first = 0
        last = 10

        cls.tree_data.Draw('tr1_cluster_energy[1]>>h_data_tr1({}, {}, {})'.format(n_bins, first, last))
        cls.tree_mc.Draw('tr1_cluster_energy[1]>>h_mc_tr1({}, {}, {})'.format(n_bins, first, last))
        h_data_tr1 = gROOT.FindObject('h_data_tr1')
        h_data_tr1.Scale(1. / h_data_tr1.GetEntries())
        h_mc_tr1 = gROOT.FindObject('h_mc_tr1')
        h_mc_tr1.Scale(1. / h_mc_tr1.GetEntries())

        cls.tree_data.Draw('tr2_cluster_energy[1]>>h_data_tr2({}, {}, {})'.format(n_bins, first, last))
        cls.tree_mc.Draw('tr2_cluster_energy[1]>>h_mc_tr2({}, {}, {})'.format(n_bins, first, last))
        h_data_tr2 = gROOT.FindObject('h_data_tr2')
        h_data_tr2.Scale(1. / h_data_tr2.GetEntries())

        h_mc_tr2 = gROOT.FindObject('h_mc_tr2')
        h_mc_tr2.Scale(1. / h_mc_tr2.GetEntries())

        h_ratio_data = h_data_tr2.Clone()
        h_ratio_data.Divide(h_data_tr1)
        h_ratio_mc = h_mc_tr2.Clone()
        h_ratio_mc.Divide(h_mc_tr1)

        h_ratio_mc.Draw('HISTOPE')
        h_ratio_data.Draw('HISTOPEsame')

        h_ratio_mc.SetMaximum(max(h_ratio_mc.GetBinContent(h_ratio_mc.GetMaximumBin()), h_ratio_data.GetBinContent(h_ratio_data.GetMaximumBin())) + 0.01)

        h_ratio_mc.SetMarkerColor(2)
        h_ratio_mc.SetMarkerStyle(20)
        h_ratio_mc.SetMarkerSize(1.5)
        h_ratio_mc.SetTitle('MC;energy, [MIP];#frac{E_{tr2}}{E_{tr1}}')

        h_ratio_data.SetMarkerStyle(20)
        h_ratio_data.SetMarkerColor(1)
        h_ratio_data.SetMarkerSize(1.5)
        h_ratio_data.SetTitle('Data')

        canvas.BuildLegend()
        h_ratio_mc.SetTitle('Trackers energy ratio')
        gStyle.SetOptStat(0)
        canvas.Write('e_ratio_tr2_tr1')

    ###

    @classmethod
    def residuals_111(cls):
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

        for event in cls.tree_data:
            if not (event.tr1_n_clusters == 1
                    and event.tr2_n_clusters == 1
                    and event.cal_n_clusters == 1
                    and 153.1 < event.cal_cluster_rho[0] < 172.3):
                continue

            rho = np.array([event.tr1_cluster_y[0], event.tr2_cluster_y[0], event.cal_cluster_y[0]])
            track = TGraphErrors(3, z, rho, z_err, rho_err)
            track.Fit('pol0', "Q")
            fit_func = track.GetFunction('pol0')

            tr1_residual = event.tr1_cluster_y[0] - fit_func.Eval(z_tr1)
            tr2_residual = event.tr2_cluster_y[0] - fit_func.Eval(z_tr2)
            cal_residual = event.cal_cluster_y[0] - fit_func.Eval(z_cal)
            h_tr1.Fill(tr1_residual)
            h_tr2.Fill(tr2_residual)
            h_cal.Fill(cal_residual)

        # for event in cls.tree_mc:
        #     if not (event.tr1_n_clusters == 1
        #             and event.tr2_n_clusters == 1
        #             and event.cal_n_clusters == 1
        #             and 153.1 < event.cal_cluster_rho[0] < 172.3):
        #         continue

        #     rho = np.array([event.tr1_cluster_y[0], event.tr2_cluster_y[0], event.cal_cluster_y[0]])
        #     track = TGraphErrors(3, z, rho, z_err, rho_err)
        #     track.Fit('pol0', "Q")
        #     fit_func = track.GetFunction('pol0')

        #     tr1_residual = event.tr1_cluster_rho[0] - fit_func.Eval(z_tr1)
        #     tr2_residual = event.tr2_cluster_rho[0] - fit_func.Eval(z_tr2)
        #     cal_residual = event.cal_cluster_rho[0] - fit_func.Eval(z_cal)
        #     h_tr1_mc.Fill(tr1_residual)
        #     h_tr2_mc.Fill(tr2_residual)
        #     h_cal_mc.Fill(cal_residual)

        # h_tr1_mc.Draw('histo')
        # h_tr1_mc.Scale(1. / h_tr1_mc.GetEntries())

        # h_tr1_mc.SetFillColor(5)

        h_tr1.Draw('histo')
        print('Tr1_event_per_event:', h_tr1.GetMean())
        h_tr1.Scale(1. / h_tr1.GetEntries())
        h_tr1.SetLineWidth(3)

        # h_tr1_mc.SetMaximum(max(h_tr1_mc.GetBinContent(h_tr1_mc.GetMaximumBin()), h_tr1.GetBinContent(h_tr1.GetMaximumBin())) + 0.01)
        canvas.BuildLegend()
        canvas.Write('tr1_res111')

        # # h_tr2_mc.Draw('histo')
        # h_tr2_mc.Scale(1. / h_tr2_mc.GetEntries())

        # h_tr2_mc.SetFillColor(5)

        h_tr2.Draw('histo')
        print('Tr2_event_per_event:', h_tr2.GetMean())

        h_tr2.Scale(1. / h_tr2.GetEntries())
        h_tr2.SetLineWidth(3)

        # h_tr2_mc.SetMaximum(max(h_tr2_mc.GetBinContent(h_tr2_mc.GetMaximumBin()), h_tr2.GetBinContent(h_tr2.GetMaximumBin())) + 0.01)
        canvas.BuildLegend()
        canvas.Write('tr2_res111')

        # h_cal_mc.Draw('histo')
        # h_cal_mc.Scale(1. / h_cal_mc.GetEntries())

        # h_cal_mc.SetFillColor(5)

        h_cal.Draw('histo')
        print('Cal_event_per_event:', h_cal.GetMean())

        h_cal.Scale(1. / h_cal.GetEntries())
        h_cal.SetLineWidth(3)

        # h_cal_mc.SetMaximum(max(h_cal_mc.GetBinContent(h_cal_mc.GetMaximumBin()), h_cal.GetBinContent(h_cal.GetMaximumBin())) + 0.01)
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
        h_mc.SetMaximum(max(h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin())) + 0.01)

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


class Calorimeter(Detector):

    def y(self):
        canvas = TCanvas('shower_y', 'shower_y', 1024, 768)
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

        h_mc.SetMaximum(max([h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin()), h_true.GetBinContent(h_true.GetMaximumBin())]) + 0.01)

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
        canvas.Write('shower_y')

    def y_loop(self):
        canvas = TCanvas('shower_y', 'shower_y', 1024, 768)
        n_bins = 180
        first = 110
        last = 200
        histo = TH1F("name", "title", 180, 110, 200)
        for event in self.tree_data:
            cluster_y = event.cal_cluster_y
            for y in cluster_y:
                histo.Fill(y)

        histo.Draw('histo')
        histo.Scale(1. / histo.GetEntries())
        histo.SetLineWidth(3)
        histo.SetTitle('Data')

        canvas.BuildLegend()
        gStyle.SetOptStat(1110)
        canvas.Write('shower_y')

    def x(self):
        canvas = TCanvas('shower_x', 'title', 1024, 768)
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

        h_mc.SetMaximum(max([h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin()), h_true.GetBinContent(h_true.GetMaximumBin())]) + 0.01)

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
        canvas.Write('shower_x')

    def energy(self):
        canvas = TCanvas('shower_energy', 'shower_energy', 1024, 768)
        n_bins = 100
        first = 0
        last = 800

        self.tree_data.Draw('cal_cluster_energy[0]>>h_data({}, {}, {})'.format(n_bins, first, last))
        self.tree_mc.Draw('cal_cluster_energy[0]>>h_mc({}, {}, {})'.format(n_bins, first, last))
        h_data = gROOT.FindObject('h_data')
        h_mc = gROOT.FindObject('h_mc')

        h_mc.Draw('histo')
        h_mc.Scale(1. / h_mc.GetEntries())
        h_data.Draw('histosame')
        h_data.Scale(1. / h_data.GetEntries())

        h_mc.SetMaximum(max(h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin())) + 0.01)

        h_mc.SetFillColor(5)
        h_mc.SetTitle('MC;energy, [MIP];#frac{N_{ev}}{N_{tot}}')

        h_data.SetLineWidth(3)
        h_data.SetTitle('Data')

        canvas.BuildLegend()
        h_mc.SetTitle('Shower energy')
        gStyle.SetOptStat(1110)
        canvas.Write('shower_energy')

    def n_pads(self):
        canvas = TCanvas('shower_n_pads', 'shower_n_pads', 1024, 768)
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

        h_mc.SetMaximum(max(h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin())) + 0.01)

        h_mc.SetFillColor(5)
        h_mc.SetTitle('MC;N_{pads};#frac{N_{ev}}{N_{tot}}')

        h_data.SetLineWidth(3)
        h_data.SetTitle('Data')

        canvas.BuildLegend()
        h_mc.SetTitle('Number of hits in shower')
        gStyle.SetOptStat(1110)
        canvas.Write('shower_n_pads')

    def n_clusters(self):
        canvas = TCanvas('shower_n_clusters', 'shower_n_clusters', 1024, 768)
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

        h_mc.SetMaximum(max(h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin())) + 0.01)

        h_mc.SetFillColor(5)
        h_mc.SetTitle('MC;N_{clusters};#frac{N_{ev}}{N_{tot}}')

        h_data.SetLineWidth(3)
        h_data.SetTitle('Data')

        canvas.BuildLegend()
        h_mc.SetTitle('Number of clusters')
        gStyle.SetOptStat(1110)
        canvas.Write('shower_n_clusters')

    def layer(self):
        canvas = TCanvas('shower_layer', 'shower_layer', 1024, 768)
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

        h_mc.SetMaximum(max(h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin())) + 0.01)

        h_mc.SetFillColor(5)
        h_mc.SetTitle('MC;layer;#frac{N_{ev}}{N_{tot}}')

        h_data.SetLineWidth(3)
        h_data.SetTitle('Data')

        canvas.BuildLegend()
        h_mc.SetTitle('Shower layer')
        gStyle.SetOptStat(1110)
        canvas.Write('shower_layer')

    def clusters_distance_ratio(self):
        canvas = TCanvas('clst1_clst2_distance_eratio', 'clst1_clst2_distance_eratio', 1024, 768)
        canvas.cd()
        gStyle.SetOptStat(1110)

        self.tree_data.Draw('nomerge.cal_cluster_energy[1]/nomerge.cal_cluster_energy[0]:(nomerge.cal_cluster_y[1]-nomerge.cal_cluster_y[0])\
                              >>h_d_to_ratio(800, -80, 80, 100, 0, 1)')

        h_d_to_ratio = gROOT.FindObject('h_d_to_ratio')
        h_d_to_ratio.SetTitle('Data, distance between clusters 1 and 2 to their energy ratio;d, [mm];Energy ratio')
        h_d_to_ratio.Draw('colz')
        canvas.Write('h_d_to_ratio')


class Tracker(Detector):

    def __init__(self, tr_idx):
        self.tr_idx = tr_idx
        super().__init__()

    def energy(self):
        name = 'tr{}_energy'.format(self.tr_idx)
        canvas = TCanvas(name, name, 1024, 768)
        langaus_fit = TF1('langaus_fit', langaufun, 0.25, 8., 4)
        langaus_fit.SetNpx(500)
        langaus_fit.SetLineColor(2)
        langaus_fit.SetMarkerStyle(1)
        langaus_fit.SetTitle('Landau*Gaus fit')
        langaus_fit.SetParNames("Width", "MP", "Area", "GSigma")
        langaus_fit.SetParameters(0.07, 1, 1., 0.13)

        self.tree_data.Draw('tr{}_cluster_energy[0]>>h_data({100, 0, 10)'.format(self.tr_idx))
        self.tree_mc.Draw('tr{}_cluster_energy[0]>>h_mc(100, 0, 10)'.format(self.tr_idx))
        h_data = gROOT.FindObject('h_data')
        h_mc = gROOT.FindObject('h_mc')

        h_mc.Draw('histo')
        h_mc.Scale(1. / h_mc.GetEntries())
        h_data.Draw('histosame')
        h_data.Scale(1. / h_data.GetEntries())

        h_data.Fit('langaus_fit')
        langaus_fit.Draw('same')

        h_mc.SetMaximum(max(h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin())) + 0.01)

        h_mc.SetFillColor(5)
        h_mc.SetTitle('MC;energy, [MIP];#frac{N_{ev}}{N_{tot}}')

        h_data.SetLineWidth(3)
        h_data.SetTitle('Data')

        canvas.BuildLegend()
        h_mc.SetTitle('Tracker{}'.format(self.tr_idx))
        gStyle.SetOptStat(1110)
        canvas.Write(name)

    def bs_energy(self):
        name = 'tr{}_bs_energy'.format(self.tr_idx)
        canvas = TCanvas(name, name, 1024, 768)
        langaus_fit = TF1('langaus_fit', langaufun, 0.5, 8., 4)
        langaus_fit.SetNpx(500)
        langaus_fit.SetLineColor(2)
        langaus_fit.SetMarkerStyle(1)
        langaus_fit.SetTitle('Landau*Gaus fit')
        langaus_fit.SetParNames("Width", "MP", "Area", "GSigma")
        langaus_fit.SetParameters(0.07, 1, 1., 0.13)

        self.tree_data.Draw('tr{}_cluster_energy[1]>>h_data(100, 0, 10)'.format(self.tr_idx))
        self.tree_mc.Draw('tr{}_cluster_energy[1]>>h_mc(100, 0, 10)'.format(self.tr_idx))
        h_data = gROOT.FindObject('h_data')
        h_mc = gROOT.FindObject('h_mc')

        h_mc.Draw('histo')
        h_mc.Scale(1. / h_mc.GetEntries())
        h_data.Draw('histosame')
        h_data.Scale(1. / h_data.GetEntries())

        h_data.Fit('langaus_fit')
        langaus_fit.Draw('same')

        h_mc.SetMaximum(max(h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin())) + 0.01)

        h_mc.SetFillColor(5)
        h_mc.SetTitle('MC;energy, [MIP];#frac{N_{ev}}{N_{tot}}')

        h_data.SetLineWidth(3)
        h_data.SetTitle('Data')

        canvas.BuildLegend()
        h_mc.SetTitle('Tracker{}'.format(self.tr_idx))
        gStyle.SetOptStat(1110)
        canvas.Write(name)

    def n_clusters(self):
        name = 'tr{}_n_clusters'.format(self.tr_idx)
        canvas = TCanvas(name, name, 1024, 768)

        self.tree_data.Draw('tr{}_n_clusters>>h_data(5, 0, 5)'.format(self.tr_idx))
        self.tree_mc.Draw('tr{}_n_clusters>>h_mc(5, 0, 5)'.format(self.tr_idx))
        h_data = gROOT.FindObject('h_data')
        h_mc = gROOT.FindObject('h_mc')

        h_mc.Draw('histo')
        h_mc.Scale(1. / h_mc.GetEntries())
        h_data.Draw('histosame')
        h_data.Scale(1. / h_data.GetEntries())

        h_mc.SetMaximum(max(h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin())) + 0.01)

        h_mc.SetFillColor(5)
        h_mc.SetTitle('MC;N_{clusters};#frac{N_{ev}}{N_{tot}}')

        h_data.SetLineWidth(3)
        h_data.SetTitle('Data')

        canvas.BuildLegend()
        h_mc.SetTitle('Tracker{}'.format(self.tr_idx))
        gStyle.SetOptStat(1110)
        canvas.Write(name)

    def efficiency(self):
        name = 'tr{}_efficiency'.format(self.tr_idx)
        canvas = TCanvas(name, name, 1024, 768)

        general_cuts = 'cal_n_clusters == 1 && tr{0}_n_clusters == 1 && cal_cluster_y[0] > 158 && cal_cluster_y[0] < 168\
                && abs(tr{0}_cluster_y - cal_cluster_y[0]) < 1.6'.format(1 if self.tr_idx == 2 else 2)

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
            selection_cuts = 'cal_n_clusters == 1 && tr{0}_n_clusters == 1 && cal_cluster_y[0] > 158 && cal_cluster_y[0] < 168\
                    && abs(tr{}_cluster_y - cal_cluster_y[0]) < 1.6 \
                    && abs(tr{1}_cluster_y - cal_cluster_y[0]) < {2}'.format(1 if self.tr_idx == 2 else 2, self.tr_idx, x)

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
        efficiency_data_graph.SetTitle('Tracker{};'.format(self.tr_idx) + '|y_{tr2} - y_{shower}|;Efficiency, %')
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

        canvas.Write(name)

    def all_clst_shower_distance(self):
        name = 'tr{}_all_clst_shower_distance'.format(self.tr_idx)
        canvas = TCanvas(name, name, 1024, 768)
        canvas.cd()

        cuts = 'cal_n_clusters == 1'
        distance = 'tr{}_cluster_y - cal_cluster_y[0]'.format(self.tr_idx)

        fit_peak = TF1("fit_peak", "gaus", -50, 30)
        fit_peak.SetMarkerStyle(1)
        fit_halo = TF1("fit_halo", "(0.5 * [0] * [1] / pi) / ((x[0] - [2]) * (x[0] - [2]) + .25 * [1] * [1])", -50, 30, 3)
        fit_halo.SetMarkerStyle(1)
        fit_total = TF1("fit_total", "gaus(0)+(0.5 * [3] * [4] / pi) / ((x[0] - [5]) * (x[0] - [5]) + .25 * [4] * [4])", -50, 30, 6)
        fit_total.SetMarkerStyle(1)

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

        self.tree_data.Draw(distance + '>>h_data(100, -70, 50)', cuts)
        self.tree_mc.Draw(distance + '>>h_mc(100, -70, 50)', cuts)

        h_data = gROOT.FindObject('h_data')
        h_mc = gROOT.FindObject('h_mc')

        h_mc.Scale(1. / h_mc.GetEntries())
        h_data.Scale(1. / h_data.GetEntries())

        h_data.Fit(fit_total)

        par123 = array.array('d', 3 * [0.])
        par456 = array.array('d', 3 * [0.])
        pars = fit_total.GetParameters()

        par123[0], par123[1], par123[2] = pars[0], pars[1], pars[2]
        par456[0], par456[1], par456[2] = pars[3], pars[4], pars[5]

        fit_peak.SetParameters(par123)
        fit_halo.SetParameters(par456)

        h_mc.Draw('histo')
        h_data.Draw('histosame')
        fit_total.Draw('same')
        # fit_peak.Draw('same')
        fit_halo.Draw('same')

        h_mc.SetMaximum(max(h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin())) + 0.01)
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
        h_mc.SetTitle('Tracker{}'.format(self.tr_idx))

        canvas.Write(name)


def main():

    cal = Calorimeter()
    cal.y_loop()
    # Detector().check_layer2_energy()
    # Detector().e_ratio_tr2_tr1()

    # cal = Calorimeter()
    # cal.energy()

    input('Yaay I am finished :3')


gROOT.SetBatch(1)
gROOT.SetStyle('ATLAS')

main()
