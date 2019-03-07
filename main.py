from ROOT import TFile, gROOT, TGraphErrors, TH1F, TGraph, TH2F, TCanvas, TF1
import numpy as np
import array


class Data():
    def __init__(self):
        self.file_data = TFile.Open('./trees/extracted_data_merge_on.root')
        self.tree_data = self.file_data.data
        self.tree_data.AddFriend('nomerge = data', './trees/extracted_data_merge_off.root')

        self.file_mc = TFile.Open('./trees/extracted_mc_merge_on.root')
        self.tree_mc = self.file_mc.mc
        self.tree_mc.AddFriend('nomerge = mc', './trees/extracted_mc_merge_on.root')

        self.output_file = TFile('RENAME.root', 'update')
        self.output_file.cd()
        self.n_events_data = self.tree_data.GetEntries()
        self.n_events_mc = self.tree_mc.GetEntries()

    def shower_distance(self, tr_number, cluster_number, tr1_n_clusters, tr2_n_clusters):
        canvas = TCanvas('canvas', 'title', 1800, 900)

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

        canvas = TCanvas('residuals_111', 'title', 1800, 1800)

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

        canvas.BuildLegend()
        canvas.Write('tr1_res111')

        h_tr2_mc.Draw('histo')
        h_tr2_mc.Scale(1. / self.n_events_mc)

        h_tr2_mc.SetFillColor(5)

        h_tr2.Draw('histosame')
        h_tr2.Scale(1. / self.n_events_data)
        h_tr2.SetLineWidth(3)

        canvas.BuildLegend()
        canvas.Write('tr2_res111')

        h_cal_mc.Draw('histo')
        h_cal_mc.Scale(1. / self.n_events_mc)

        h_cal_mc.SetFillColor(5)

        h_cal.Draw('histosame')
        h_cal.Scale(1. / self.n_events_data)
        h_cal.SetLineWidth(3)

        canvas.BuildLegend()
        canvas.Write('cal_res111')

    def beam_pos_cal(self):
        cut_option = ['no_cut', 'cuts']
        cuts = ['', 'cal_n_clusters == 1 && cal_cluster_energy[0] > 100 && cal_cluster_rho[0] > 153.1 && cal_cluster_rho[0] < 172.3']
        for idx, cut in enumerate(cuts):
            canvas = TCanvas('beam_pos_cal_' + cut_option[idx], 'title', 1800, 900)
            n_bins = 200
            first = 120
            last = 195.2

            self.tree_data.Draw('cal_cluster_rho[0]>>h_data({}, {}, {})'.format(n_bins, first, last), cut)
            self.tree_mc.Draw('cal_cluster_rho[0]>>h_mc({}, {}, {})'.format(n_bins, first, last), cut)
            self.tree_mc.Draw('cal_true_hit_rho>>h_true({}, {}, {})'.format(n_bins, first, last), cut)
            h_mc = gROOT.FindObject('h_mc')
            h_data = gROOT.FindObject('h_data')
            h_true = gROOT.FindObject('h_true')

            h_mc.Draw('histo')
            h_mc.Scale(1. / self.n_events_mc)
            h_data.Draw('histosame')
            h_data.Scale(1. / self.n_events_data)
            h_true.Draw('histosame')
            h_true.Scale(1. / self.n_events_mc)
            h_mc.SetMaximum(max([h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin()), h_true.GetBinContent(h_true.GetMaximumBin())]))

            h_mc.SetFillColor(5)
            h_mc.SetTitle('MC')

            h_data.SetLineWidth(3)
            h_data.SetTitle('Data')

            h_true.SetLineWidth(3)
            h_true.SetLineColor(2)
            h_true.SetTitle('Generated')

            canvas.BuildLegend()
            h_mc.SetTitle('beam_profile_' + cut_option[idx])

            canvas.Write('beam_profile_' + cut_option[idx])

    def distance_between_clusters(self, det, clst1, clst2):
        canvas = TCanvas('canvas', 'distance between clusters', 1800, 1800)
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

    def all_events_shower_distance(self):
        canvas = TCanvas('canvas', 'title', 1800, 900)

        h_name = 'all_events_tr2_clst_distance_to_shower'

        cuts = 'cal_n_clusters == 1 && cal_cluster_rho[0] > 153.1 && cal_cluster_rho[0] < 172.3 && cal_cluster_energy[0] > 100'
        distance = 'tr2_cluster_rho - cal_cluster_rho[0]'

        fit_total = TF1('fit_total', 'gaus(0) + gaus(3)', -60, 40)
        fit_peak = TF1('fit_peak', 'gaus', -60, 40)
        fit_halo = TF1('fit_halo', 'gaus', -60, 40)
        fit_total.SetParLimits(2, 0.4, 1.)
        fit_total.SetParLimits(5, 12., 30.)
        fit_total.SetNpx(500)
        fit_peak.SetNpx(500)
        fit_halo.SetNpx(500)

        self.tree_data.Draw(distance + '>>' + h_name + '(300, -60, 40)', cuts)
        self.tree_mc.Draw(distance + '>>' + h_name + '_mc' + '(300, -60, 40)', cuts)

        h_data = gROOT.FindObject(h_name)
        h_mc = gROOT.FindObject(h_name + '_mc')

        h_mc.Scale(1. / self.n_events_mc)
        h_data.Scale(1. / self.n_events_data)

        h_mc.Fit(fit_total)

        par1 = array.array('d', 3 * [0.])
        par2 = array.array('d', 3 * [0.])
        pars = fit_total.GetParameters()

        par1[0], par1[1], par1[2] = pars[0], pars[1], pars[2]
        par2[0], par2[1], par2[2] = pars[3], pars[4], pars[5]

        fit_peak.SetParameters(par1)
        fit_halo.SetParameters(par2)

        h_mc.Draw('histo')
        h_data.Draw('histosame')
        fit_total.Draw('same')
        fit_peak.Draw('same')
        #fit_halo.Draw('same')

        h_mc.SetMaximum(max(h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin())))
        h_mc.SetTitle(h_name + '_mc')
        h_mc.GetXaxis().SetTitle('Distance to the shower, mm')
        h_mc.GetYaxis().SetTitle('N_{events} normalized')
        h_mc.SetFillColor(5)

        h_data.SetTitle(h_name)
        h_data.SetLineWidth(3)

        fit_total.SetLineColor(2)
        fit_total.SetTitle('Total fit')
        fit_peak.SetLineColor(6)
        fit_peak.SetTitle('Peak fit')
        fit_halo.SetLineColor(7)
        fit_halo.SetTitle('Halo fit')

        canvas.SetLogy()
        canvas.BuildLegend()

        canvas.Write(h_name)

    def n_clusters(self, det):
        canvas = TCanvas('canvas', 'title', 1800, 900)
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
            canvas = TCanvas('canvas', 'title', 1800, 900)
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
        canvas = TCanvas('merge_pos_change', 'title', 1800, 900)

        self.tree_data.Draw('(nomerge.cal_cluster_rho[0] - cal_cluster_rho[0])>>h_data')
        self.tree_mc.Draw('(nomerge.cal_cluster_y[0] - cal_cluster_y[0])>>h_mc')

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

        canvas.Write('merge_pos_shift')

    def backscattered_tracks(self):
        # z_cal_entrace = 23 * 4.5
        cuts = 'cal_n_clusters == 1 && cal_cluster_rho[0] > 153.1 && cal_cluster_rho[0] < 172.3 && cal_cluster_energy[0] > 100'

        self.tree_data.Draw('(tr1_cluster_y[1]+23*4.5*(tr2_cluster_y[1]-tr1_cluster_y[1])/22.5):\
                             (tr1_cluster_x[1]+23*4.5*(tr2_cluster_x[1]-tr1_cluster_x[1])/22.5)\
                              >>h_tracks(100, -51, 51, 100, 0, 196)', cuts)

        h_tracks = gROOT.FindObject('h_tracks')
        h_tracks.SetTitle('Data, backscattered_tracks; x, [mm]; y, [mm]')

        h_tracks.Write()

        self.tree_mc.Draw('(tr1_cluster_y[1]+23*4.5*(tr2_cluster_y[1]-tr1_cluster_y[1])/22.5):\
                             (tr1_cluster_x[1]+23*4.5*(tr2_cluster_x[1]-tr1_cluster_x[1])/22.5)\
                              >>h_tracks_mc(100, -51, 51, 100, 0, 196)', cuts)

        h_tracks_mc = gROOT.FindObject('h_tracks_mc')
        h_tracks_mc.SetTitle('MC, backscattered_tracks; x, [mm]; y, [mm]')

        h_tracks_mc.Write()


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


def tr2_residual_efficiency(tree, output_file):
    x_arr = np.arange(0, 15, 0.2)
    eff = np.zeros_like(x_arr)
    eff_tot = 0

    tr1_z = 0.
    tr2_z = 5 * 4.5
    cal_z = 25 * 4.5
    fit_x = np.array([tr1_z, cal_z])
    fit_x_err = np.array([4.5 / 2, 13.5])
    fit_y_err = np.array([0.9, 0.9])
    for event in tree:
        if (event.tr1_n_clusters >= 1 and event.cal_n_clusters == 1
           and 153.1 < event.cal_cluster_rho[0] < 172.3
           and -4.4 < event.tr1_cluster_rho[0] - event.cal_cluster_rho[0] < 4.4):
            eff_tot += 1
            fit_y = np.array([event.tr1_cluster_rho[0], event.cal_cluster_rho[0]])
            track = TGraphErrors(2, fit_x, fit_y, fit_x_err, fit_y_err)
            track.Fit('pol0', "Q")
            fit_func = track.GetFunction('pol0')
            for idx, x in enumerate(x_arr):
                    if event.tr2_n_clusters >= 1 and -x < event.tr2_cluster_rho[0] - fit_func.Eval(tr2_z) < x:
                        eff[idx] += 1

    eff_graph = TGraph(len(x_arr), x_arr, eff / eff_tot)
    eff_graph.SetName('tr2_eff_res')
    eff_graph.SetTitle('Tracker2 Efficiency, residual; d_{tr2-fit}; Eff, %')
    eff_graph.SetLineWidth(2)
    output_file.cd()
    eff_graph.Write()


def tr2_distance_efficiency(tree, output_file):
    x_arr = np.arange(0, 15, 0.2)
    eff = np.zeros_like(x_arr)
    eff_tot = 0

    for event in tree:
        if (event.tr1_n_clusters >= 1 and event.cal_n_clusters == 1
           and 153.1 < event.cal_cluster_rho[0] < 172.3
           and -4.4 < event.tr1_cluster_rho[0] - event.cal_cluster_rho[0] < 4.4):
            eff_tot += 1
            for idx, x in enumerate(x_arr):
                    if event.tr2_n_clusters >= 1 and -x < event.tr2_cluster_rho[0] - event.cal_cluster_rho[0] < x:
                        eff[idx] += 1

    eff_graph = TGraph(len(x_arr), x_arr, eff / eff_tot)
    eff_graph.SetName('tr2_eff_dist')
    eff_graph.SetTitle('Tracker2 Efficiency, distance; d_{tr2-shower}; Eff, %')
    eff_graph.SetLineWidth(2)
    output_file.cd()
    eff_graph.Write()


def tr1_residual_efficiency(tree, output_file):
    x_arr = np.arange(0, 15, 0.2)
    eff = np.zeros_like(x_arr)
    eff_tot = 0

    tr1_z = 0.
    tr2_z = 5 * 4.5
    cal_z = 25 * 4.5
    fit_x = np.array([tr2_z, cal_z])
    fit_x_err = np.array([4.5 / 2, 13.5])
    fit_y_err = np.array([0.9, 0.9])
    for event in tree:
        if (event.tr2_n_clusters >= 1 and event.cal_n_clusters == 1
           and 153.1 < event.cal_cluster_rho[0] < 172.3
           and -4.4 < event.tr2_cluster_rho[0] - event.cal_cluster_rho[0] < 4.4):
            eff_tot += 1
            fit_y = np.array([event.tr2_cluster_rho[0], event.cal_cluster_rho[0]])
            track = TGraphErrors(2, fit_x, fit_y, fit_x_err, fit_y_err)
            track.Fit('pol0', "Q")
            fit_func = track.GetFunction('pol0')
            for idx, x in enumerate(x_arr):
                    if event.tr1_n_clusters >= 1 and -x < event.tr1_cluster_rho[0] - fit_func.Eval(tr1_z) < x:
                        eff[idx] += 1

    eff_graph = TGraph(len(x_arr), x_arr, eff / eff_tot)
    eff_graph.SetName('tr1_eff_res')
    eff_graph.SetTitle('Tracker1 Efficiency, residual; d_{tr1-fit}; Eff, %')
    eff_graph.SetLineWidth(2)
    output_file.cd()
    eff_graph.Write()


def tr1_distance_efficiency(tree, output_file):
    x_arr = np.arange(0, 15, 0.2)
    eff = np.zeros_like(x_arr)
    eff_tot = 0

    for event in tree:
        if (event.tr2_n_clusters >= 1 and event.cal_n_clusters == 1
           and 153.1 < event.cal_cluster_rho[0] < 172.3
           and -4.4 < event.tr2_cluster_rho[0] - event.cal_cluster_rho[0] < 4.4):
            eff_tot += 1
            for idx, x in enumerate(x_arr):
                    if event.tr1_n_clusters >= 1 and -x < event.tr1_cluster_rho[0] - event.cal_cluster_rho[0] < x:
                        eff[idx] += 1

    eff_graph = TGraph(len(x_arr), x_arr, eff / eff_tot)
    eff_graph.SetName('tr1_eff_dist')
    eff_graph.SetTitle('Tracker1 Efficiency, distance; d_{tr1-shower}; Eff, %')
    eff_graph.SetLineWidth(2)
    output_file.cd()
    eff_graph.Write()


def main():

    data = Data()
    data.residuals_111()
    data.beam_pos_cal()
    data.resolution('cal')
    data.resolution('tr1')
    data.resolution('tr2')
    data.all_events_shower_distance()
    data.n_clusters('cal')
    data.clusters_energy()
    data.merge_pos_change()
    data.backscattered_tracks()

    input('Yaay I am finished :3')


gROOT.SetBatch(1)

main()
