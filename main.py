from ROOT import TFile, gROOT, TGraphErrors, TH1F, TGraph, TH2F, TCanvas, TF1
import numpy as np
import array


def residuals_111_before_alignment(tree, output_file):
    '''
    Mean of residuals before misalignment:
    Tr1: -0.1897
    Tr2: 0.9405
    Cal: -0.7501
    '''
    tr1_z = 0.
    tr2_z = 5 * 4.5
    cal_z = 25 * 4.5
    x = np.array([tr1_z, tr2_z, cal_z])
    h_tr1_res = TH1F('h_tr1_res111', 'Tracker1 residuals', 200, -20, 20)
    h_tr2_res = TH1F('h_tr2_res111', 'Tracker2 residuals', 200, -20, 20)
    h_cal_res = TH1F('h_cal_res111', 'Calorimeter residuals', 200, -20, 20)

    for event in tree:
        if (event.tr1_n_clusters == 1 and event.tr2_n_clusters == 1 and event.cal_n_clusters == 1
           and 153.1 < event.cal_cluster_pad[0] * 1.8 + 0.9 + 80. < 172.3):
            y = np.array([event.tr1_cluster_pad[0] * 1.8 + 0.9 + 80., event.tr2_cluster_pad[0] * 1.8 + 0.9 + 80., event.cal_cluster_pad[0] * 1.8 + 0.9 + 80.])
            x_err = np.array([4.5 / 2, 4.5 / 2, 13.5])
            y_err = np.array([0.9, 0.9, 0.9])
            track = TGraphErrors(3, x, y, x_err, y_err)
            track.Fit('pol0', "Q")
            fit_func = track.GetFunction('pol0')

            tr1_residual = event.tr1_cluster_pad[0] * 1.8 + 0.9 + 80. - fit_func.Eval(tr1_z)
            tr2_residual = event.tr2_cluster_pad[0] * 1.8 + 0.9 + 80. - fit_func.Eval(tr2_z)
            cal_residual = event.cal_cluster_pad[0] * 1.8 + 0.9 + 80. - fit_func.Eval(cal_z)
            h_tr1_res.Fill(tr1_residual)
            h_tr2_res.Fill(tr2_residual)
            h_cal_res.Fill(cal_residual)

    output_file.cd()
    h_tr1_res.Write()
    h_tr2_res.Write()
    h_cal_res.Write()


class Data():
    def __init__(self):
        self.file_data = TFile.Open('./trees/extracted_data_not_merged.root')
        self.tree_data = self.file_data.data
        self.tree_data.AddFriend('merged = data', './trees/extracted_data_merged.root')

        self.file_mc = TFile.Open('./trees/extracted_mc_not_merged.root')
        self.tree_mc = self.file_mc.mc
        self.tree_mc.AddFriend('merged = mc', './trees/extracted_mc_merged.root')

        self.output_file = TFile('RENAME.root', 'recreate')
        #self.output_file_merged = TFile('RENAME_merged.root', 'recreate')

        self.n_events_data = self.tree_data.GetEntries()
        self.n_events_mc = self.tree_mc.GetEntries()
        self.n_events_data_passed = 0
        self.n_events_mc_passed = 0

    def shower_distance(self, tr_number, cluster_number, tr1_n_clusters, tr2_n_clusters):
        self.output_file.cd()
        canvas = TCanvas('canvas', 'title', 1800, 900)

        h_name = 'tr{}_clst{}_shower_distance_{}{}1'.format(tr_number, cluster_number, tr1_n_clusters, tr2_n_clusters)

        distance = 'tr{}_cluster_rho[{}] - cal_cluster_rho[0]'.format(tr_number, cluster_number - 1)

        cuts = 'tr1_n_clusters == {} && tr2_n_clusters == {}'.format(tr1_n_clusters, tr2_n_clusters)

        self.tree_data.Draw(distance + '>>' + h_name + '(100, -60, 40)', cuts)
        self.tree_mc.Draw(distance + '>>' + h_name + '_mc' + '(100, -60, 40)', cuts)

        h_data = gROOT.FindObject(h_name)
        h_mc = gROOT.FindObject(h_name + '_mc')

        self.n_events_data_passed = h_data.GetEntries()
        self.n_events_mc_passed = h_mc.GetEntries()

        h_mc.Draw('histo')
        h_mc.Scale(1. / self.n_events_mc)
        h_data.Draw('histosame')
        h_data.Scale(1. / self.n_events_data)

        h_mc.SetMaximum(max(h_mc.GetBinContent(h_mc.GetMaximumBin()), h_data.GetBinContent(h_data.GetMaximumBin())))
        h_mc.SetTitle(h_name + '_mc')
        h_mc.GetXaxis().SetTitle('Distance to the shower, mm')
        h_mc.GetYaxis().SetTitle('N_{events} normalized')
        h_mc.SetFillColor(5)

        h_data.SetTitle(h_name)
        h_data.SetLineWidth(3)

        canvas.SetLogy()

        canvas.BuildLegend()

        canvas.Write(h_name)
        canvas.Print(h_name + '.png')

    def residuals_111(self):
        z_tr1 = 0.
        z_tr2 = 5 * 4.5
        z_cal = 25 * 4.5
        z = np.array([z_tr1, z_tr2, z_cal])
        z_err = np.array([4.5 / 2, 4.5 / 2, 13.5])
        rho_err = np.array([0.9, 0.9, 0.9])

        canvas = TCanvas('canvas', 'Residuals 111', 1800, 1800)

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
        self.output_file.cd()
        cut_option = ['no_cut', 'position_cut_y_from_155_to_170', 'position_cutted_events']
        cuts = ['', 'cal_cluster_y[0] < 170 && cal_cluster_y[0] > 155', 'cal_cluster_y[0] > 170 || cal_cluster_y[0] < 155']
        for idx, cut in enumerate(cuts):
            canvas = TCanvas(cut_option[idx], 'beam_pos_cal', 1800, 1800)
            n_bins = 200
            first = 80
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
            h_mc.SetTitle('beam_profile_merged_' + cut_option[idx])

            canvas.Write('beam_profile_merged_' + cut_option[idx])

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

    def cluster_minus_generated_pos(self, det):
        canvas = TCanvas('canvas', 'title', 1800, 1800)
        n_bins = 120
        first = -30
        last = 30

        self.tree_mc.Draw(det + '_cluster_rho[0] - ' + det + '_true_hit_rho>>h_mc({}, {}, {})'.format(n_bins, first, last))
        h_mc = gROOT.FindObject('h_mc')

        h_mc.Draw('histo')
        h_mc.Scale(1. / self.n_events_mc)

        h_mc.SetFillColor(5)
        h_mc.SetTitle('MC')

        canvas.BuildLegend()
        canvas.Write(det + '_cluster_to_expected_distance')

    def all_events_shower_distance(self):
        self.output_file.cd()
        canvas = TCanvas('canvas', 'title', 1800, 900)

        h_name = 'all_events_tr2_clst_distance_to_shower'

        distance = 'tr2_cluster_rho - cal_cluster_rho[0]'

        fit_total = TF1('fit_total', 'gaus(0) + gaus(3)', -60, 40)
        fit_peak = TF1('fit_peak', 'gaus', -60, 40)
        fit_halo = TF1('fit_halo', 'gaus', -60, 40)
        fit_total.SetParLimits(2, 0.4, 1.)
        fit_total.SetParLimits(5, 12., 30.)
        fit_total.SetNpx(500)
        fit_peak.SetNpx(500)
        fit_halo.SetNpx(500)

        self.tree_data.Draw(distance + '>>' + h_name + '(300, -60, 40)')
        self.tree_mc.Draw(distance + '>>' + h_name + '_mc' + '(300, -60, 40)')

        h_data = gROOT.FindObject(h_name)
        h_mc = gROOT.FindObject(h_name + '_mc')

        self.n_events_data_passed = h_data.GetEntries()
        self.n_events_mc_passed = h_mc.GetEntries()

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
        fit_halo.Draw('same')

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
        canvas.Print(h_name + '.png')

    def n_clusters(self, det):
        canvas = TCanvas('canvas', 'title', 1800, 1800)
        n_bins = 15
        first = 0
        last = 15

        self.tree_data.Draw(det + '_n_clusters>>h_data({}, {}, {})'.format(n_bins, first, last))
        self.tree_mc.Draw(det + '_n_clusters>>h_mc({}, {}, {})'.format(n_bins, first, last))

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
            canvas = TCanvas('canvas', 'title', 1800, 1800)
            n_bins = 500
            first = 0
            last = 50

            tree.Draw('cal_cluster_energy[0]*0.0885*0.917112>>h_clst1_energy({}, {}, {})'.format(n_bins, first, last))
            tree.Draw('cal_cluster_energy[1]*0.0885*0.917112>>h_clst2_energy({}, {}, {})'.format(n_bins, first, last))
            tree.Draw('cal_cluster_energy[2]*0.0885*0.917112>>h_clst3_energy({}, {}, {})'.format(n_bins, first, last))

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

            h_clst1_energy.GetXaxis().SetTitle('Energy, [MeV]')
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

    def cluster1_cluster2_energies(self):

        self.tree_data.Draw('cal_cluster_energy[1]*0.0885*0.917112:cal_cluster_energy[0]*0.0885*0.917112>>h_clsters_energy', 'cal_cluster_energy[1]*0.0885*0.917112>5')
        h_clsters_energy = gROOT.FindObject('h_clsters_energy')
        h_clsters_energy.Draw('colz')
        h_clsters_energy.SetTitle('clst1 clst2 energies;Cluster1 Energy, [MeV]; Cluster2 Energy, [MeV]')

        #h_clsters_energy.GetXaxis().SetTitle('Energy, [MeV]')
        #h_clsters_energy.GetYaxis().SetTitle('N_{events} normalized')
        #h_clsters_energy.SetLineColor(1)
        #h_clsters_energy.SetLineWidth(3)

        #canvas.BuildLegend()
        h_clsters_energy.Write('clusters12_energies_data')

    def merge_pos_change(self):
        self.output_file.cd()
        canvas = TCanvas('canvas', 'beam_pos_cal', 1800, 1800)

        self.tree_data.Draw('(merged.cal_cluster_y[0] - cal_cluster_y[0]):cal_cluster_energy[1]*0.0885*0.917112>>h_data')
        h_data = gROOT.FindObject('h_data')
        h_data.SetTitle('Position shift vs cluster2 energy; cluster2 energy, [MeV];position shift, [mm]')
        h_data.Write('data_position_shiftvsEclst2')

        self.tree_mc.Draw('(merged.cal_cluster_y[0] - cal_cluster_y[0]):cal_cluster_energy[1]*0.0885*0.917112>>h_mc')
        h_mc = gROOT.FindObject('h_mc')
        h_mc.SetTitle('Position shift vs cluster2 energy; cluster2 energy, [MeV];position shift, [mm]')
        h_mc.Write('MC_position_shiftvsEclst2')


def shower_distances(data):
    '''
    Calculate efficiency for Tr2 if Tr1 has signal
    in range of +-"sigma" = 4.4 mm within shower in cal
    '''

    # tr_number, cluster_number, tr1_n_clusters, tr2_n_clusters
    # 111
    data.shower_distance(1, 1, 1, 1)
    print('Data Events 111:', data.n_events_data_passed * 100 / data.n_events_data)
    print('MC Events 111:', data.n_events_mc_passed * 100 / data.n_events_mc)
    data.shower_distance(2, 1, 1, 1)
    # 121
    data.shower_distance(1, 1, 1, 2)
    print('Data Events 121:', data.n_events_data_passed * 100 / data.n_events_data)
    print('MC Events 121:', data.n_events_mc_passed * 100 / data.n_events_mc)
    data.shower_distance(2, 1, 1, 2)
    data.shower_distance(2, 2, 1, 2)
    # 211
    data.shower_distance(1, 1, 2, 1)
    print('Data Events 211:', data.n_events_data_passed * 100 / data.n_events_data)
    print('MC Events 211:', data.n_events_mc_passed * 100 / data.n_events_mc)
    data.shower_distance(1, 2, 2, 1)
    data.shower_distance(2, 1, 2, 1)
    # 221
    data.shower_distance(1, 1, 2, 2)
    print('Data Events 221:', data.n_events_data_passed * 100 / data.n_events_data)
    print('MC Events 221:', data.n_events_mc_passed * 100 / data.n_events_mc)
    data.shower_distance(1, 2, 2, 2)
    data.shower_distance(2, 1, 2, 2)
    data.shower_distance(2, 2, 2, 2)


def shower_distances_energy(tree, output_file):
    '''
    Calculate efficiency for Tr2 if Tr1 has signal
    in range of +-"sigma" = 4.4 mm within shower in cal
    '''
    output_file.cd()

    # ##################################### 111
    tree.Draw('tr1_cluster_energy[0]:(tr1_cluster_rho[0] - cal_cluster_rho[0])>>h_tr1_shower_distance_111',
              'cal_n_clusters == 1 && 153.1 < cal_cluster_rho[0] && cal_cluster_rho[0] < 172.3 \
              && tr1_n_clusters == 1 && tr2_n_clusters == 1')

    histo = gROOT.FindObject('h_tr1_shower_distance_111')
    histo.SetTitle('Tracker1-Shower distance_111; distance to shower, mm; Cluster energy')
    histo.Write()

    tree.Draw('tr2_cluster_energy[0]:(tr2_cluster_rho[0] - cal_cluster_rho[0])>>h_tr2_shower_distance_111',
              'cal_n_clusters == 1 && 153.1 < cal_cluster_rho[0] && cal_cluster_rho[0] < 172.3 \
              && tr1_n_clusters == 1 && tr2_n_clusters == 1')

    histo = gROOT.FindObject('h_tr2_shower_distance_111')
    histo.SetTitle('Tracker2-Shower distance_111; distance to shower, mm; Cluster energy')
    histo.Write()

    # ##################################### 121
    tree.Draw('tr1_cluster_energy[0]:(tr1_cluster_rho[0] - cal_cluster_rho[0])>>h_tr1_shower_distance_121',
              'cal_n_clusters == 1 && 153.1 < cal_cluster_rho[0] && cal_cluster_rho[0] < 172.3 \
              && tr1_n_clusters == 1 && tr2_n_clusters == 2')

    histo = gROOT.FindObject('h_tr1_shower_distance_121')
    histo.SetTitle('Tracker1-Shower distance_121; distance to shower, mm; Cluster energy')
    histo.Write()

    tree.Draw('tr2_cluster_energy[0]:(tr2_cluster_rho[0] - cal_cluster_rho[0])>>h_tr2_clst1_shower_distance_121',
              'cal_n_clusters == 1 && 153.1 < cal_cluster_rho[0] && cal_cluster_rho[0] < 172.3 \
              && tr1_n_clusters == 1 && tr2_n_clusters == 2')

    histo = gROOT.FindObject('h_tr2_clst1_shower_distance_121')
    histo.SetTitle('Tracker2 Clst1-Shower distance_121; distance to shower, mm; Cluster energy')
    histo.Write()

    tree.Draw('tr2_cluster_energy[1]:(tr2_cluster_rho[1] - cal_cluster_rho[0])>>h_tr2_clst2_shower_distance_121',
              'cal_n_clusters == 1 && 153.1 < cal_cluster_rho[0] && cal_cluster_rho[0] < 172.3 \
              && tr1_n_clusters == 1 && tr2_n_clusters == 2')

    histo = gROOT.FindObject('h_tr2_clst2_shower_distance_121')
    histo.SetTitle('Tracker2 Clst2-Shower distance_121; distance to shower, mm; Cluster energy')
    histo.Write()
    # ##################################### 211
    tree.Draw('tr1_cluster_energy[0]:(tr1_cluster_rho[0] - cal_cluster_rho[0])>>h_tr1_clst1_shower_distance_211',
              'cal_n_clusters == 1 && 153.1 < cal_cluster_rho[0] && cal_cluster_rho[0] < 172.3 \
              && tr1_n_clusters == 2 && tr2_n_clusters == 1')

    histo = gROOT.FindObject('h_tr1_clst1_shower_distance_211')
    histo.SetTitle('Tracker1 Clst1-Shower distance_211; distance to shower, mm; Cluster energy')
    histo.Write()

    tree.Draw('tr1_cluster_energy[1]:(tr1_cluster_rho[1] - cal_cluster_rho[0])>>h_tr1_clst2_shower_distance_211',
              'cal_n_clusters == 1 && 153.1 < cal_cluster_rho[0] && cal_cluster_rho[0] < 172.3 \
              && tr1_n_clusters == 2 && tr2_n_clusters == 1')

    histo = gROOT.FindObject('h_tr1_clst2_shower_distance_211')
    histo.SetTitle('Tracker1 Clst2-Shower distance_211; distance to shower, mm; Cluster energy')
    histo.Write()

    tree.Draw('tr2_cluster_energy[0]:(tr2_cluster_rho[0] - cal_cluster_rho[0])>>h_tr2_shower_distance_211',
              'cal_n_clusters == 1 && 153.1 < cal_cluster_rho[0] && cal_cluster_rho[0] < 172.3 \
              && tr1_n_clusters == 2 && tr2_n_clusters == 1')

    histo = gROOT.FindObject('h_tr2_shower_distance_211')
    histo.SetTitle('Tracker2-Shower distance_211; distance to shower, mm; Cluster energy')
    histo.Write()
    # ##################################### 221
    tree.Draw('tr1_cluster_energy[0]:(tr1_cluster_rho[0] - cal_cluster_rho[0])>>h_tr1_clst1_shower_distance_221',
              'cal_n_clusters == 1 && 153.1 < cal_cluster_rho[0] && cal_cluster_rho[0] < 172.3 \
              && tr1_n_clusters == 2 && tr2_n_clusters == 2')

    histo = gROOT.FindObject('h_tr1_clst1_shower_distance_221')
    histo.SetTitle('Tracker1 Clst1-Shower distance_221; distance to shower, mm; Cluster energy')
    histo.Write()

    tree.Draw('tr1_cluster_energy[1]:(tr1_cluster_rho[1] - cal_cluster_rho[0])>>h_tr1_clst2_shower_distance_221',
              'cal_n_clusters == 1 && 153.1 < cal_cluster_rho[0] && cal_cluster_rho[0] < 172.3 \
              && tr1_n_clusters == 2 && tr2_n_clusters == 2')

    histo = gROOT.FindObject('h_tr1_clst2_shower_distance_221')
    histo.SetTitle('Tracker1 Clst2-Shower distance_221; distance to shower, mm; Cluster energy')
    histo.Write()

    tree.Draw('tr2_cluster_energy[0]:(tr2_cluster_rho[0] - cal_cluster_rho[0])>>h_tr2_clst1_shower_distance_221',
              'cal_n_clusters == 1 && 153.1 < cal_cluster_rho[0] && cal_cluster_rho[0] < 172.3 \
              && tr1_n_clusters == 2 && tr2_n_clusters == 2')

    histo = gROOT.FindObject('h_tr2_clst1_shower_distance_221')
    histo.SetTitle('Tracker2 Clst1-Shower distance_221; distance to shower, mm; Cluster energy')
    histo.Write()

    tree.Draw('tr2_cluster_energy[1]:(tr2_cluster_rho[1] - cal_cluster_rho[0])>>h_tr2_clst2_shower_distance_221',
              'cal_n_clusters == 1 && 153.1 < cal_cluster_rho[0] && cal_cluster_rho[0] < 172.3 \
              && tr1_n_clusters == 2 && tr2_n_clusters == 2')

    histo = gROOT.FindObject('h_tr2_clst2_shower_distance_221')
    histo.SetTitle('Tracker2 Clst2-Shower distance_221; distance to shower, mm; Cluster energy')
    histo.Write()


def plot_backscattered_tracks(tree, output_file):
    h_tracks = TH2F('h_tracks', 'Backscattered tracks', 100, -51, 51, 100, 0, 196)
    z_cal_entrace = 23 * 4.5
    output_file.cd()
    for event in tree:
        if (event.tr1_n_clusters == 2 and event.tr2_n_clusters == 2 and event.cal_n_clusters == 1
           and 153.1 < event.cal_cluster_rho[0] < 172.3):
            x0 = event.tr1_cluster_x[1]
            kx = (event.tr2_cluster_x[1] - x0) / 22.5  # mm
            x = x0 + z_cal_entrace * kx

            y0 = event.tr1_cluster_y[1]
            ky = (event.tr2_cluster_y[1] - y0) / 22.5  # mm
            y = y0 + z_cal_entrace * ky

            h_tracks.Fill(x, y)
    h_tracks.Draw()
    h_tracks.Write()


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


def merge_pos_change_obsolete(input_merged_file, input_not_merged_file, output_file):
    tree_merged = input_merged_file.data
    tree_not_merged = input_not_merged_file.data

    h_pos_shift = TH2F('h_pos_shift_vs_clst2_energy', 'Merge shift position', 100, -10, 10, 100, 0, 20)

    merged_pos = np.array([])
    not_merged_pos = np.array([])
    clst2_energy = np.array([])

    print('I am here 1')
    for event1 in tree_merged:
        if event1.cal_n_clusters >= 1:
            merged_pos = np.append(merged_pos, event1.cal_cluster_y[0])

    print('I am here 2')
    for event2 in tree_not_merged:
        if event2.cal_n_clusters >= 1:
            not_merged_pos = np.append(not_merged_pos, event2.cal_cluster_y[0])
            if event1.cal_n_clusters >= 2:
                clst2_energy = np.append(clst2_energy, event1.cal_cluster_energy[1])
            else:
                clst2_energy = np.append(clst2_energy, 0)

    print('I am here 3')
    for idx, pos in enumerate(merged_pos):
        h_pos_shift.Fill(not_merged_pos[idx] - pos, clst2_energy[idx])

    output_file.cd()
    h_pos_shift.Write()


def merge_cluster_energy_change(input_merged_file, input_not_merged_file, output_file):
    tree_merged = input_merged_file.data
    tree_not_merged = input_not_merged_file.data

    h_energy_diff = TH2F('h_energy_diff', '', 200, 0, 70, 200, 0, 70)
    h_energy_diff.SetTitle('Merge vs not merged energy;Merged Energy, [MeV]; Not Merged Energy, [MeV]')

    merged_energy = np.array([])
    not_merged_energy = np.array([])

    print('I am here 1')
    for event1 in tree_merged:
        if event1.cal_n_clusters >= 1:
            merged_energy = np.append(merged_energy, event1.cal_cluster_energy[0] * 0.0885 * 0.917112)

    print('I am here 2')
    for event2 in tree_not_merged:
        if event2.cal_n_clusters >= 1:
            not_merged_energy = np.append(not_merged_energy, event2.cal_cluster_energy[0] * 0.0885 * 0.917112)

    print('I am here 3')
    for idx, energy in enumerate(merged_energy):
        h_energy_diff.Fill(not_merged_energy[idx], energy)

    output_file.cd()
    h_energy_diff.Write()


def main():

    data = Data()
    data.beam_pos_cal()
    data.cluster_minus_generated_pos('cal')
    shower_distances(data)
    data.all_events_shower_distance()
    data.merge_pos_change()

    # data.cluster1_cluster2_energies()
    # merge_cluster_energy_change(file_data_merged, file_data, output_file)

    # data.distance_between_clusters('tr1', 0, 1)

    input('Yaay I am finished :3')


gROOT.SetBatch(1)

main()

'''
TO DO LIST:

'''
