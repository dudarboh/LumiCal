
from ROOT import TFile, gROOT, gDirectory, TGraphErrors, TCanvas, gStyle, nullptr, TH2F
import numpy as np


class Detector:
    file_5gev = TFile.Open('./result_trees/extracted_run741_5gev.root', 'read')
    tree_5gev = file_5gev.data

    file_4gev = TFile.Open('./result_trees/extracted_run742_4gev.root', 'read')
    tree_4gev = file_4gev.data

    file_3gev = TFile.Open('./result_trees/extracted_run745_3gev.root', 'read')
    tree_3gev = file_3gev.data

    file_2gev = TFile.Open('./result_trees/extracted_run747_2gev.root', 'read')
    tree_2gev = file_2gev.data

    file_1gev = TFile.Open('./result_trees/extracted_run750_1gev.root', 'read')
    tree_1gev = file_1gev.data

    trees = [tree_1gev, tree_2gev, tree_3gev, tree_4gev, tree_5gev]

    output_file = TFile('histos.root', 'recreate')

    def __init__(self):
        self.output_file.cd()

    @classmethod
    def check_alignment(cls, tree_data):
        x = np.array([0., 5 * 4.5, 25 * 4.5])

        cuts = "tr1_n_clusters == 1 && tr2_n_clusters == 1 && cal_n_clusters == 1 && 153.1 < cal_cluster_y[0] && cal_cluster_y[0] < 172.3"

        tree_data.Draw('tr1_cluster_y[0]>>h_tr1(64, 80, 195.2)', cuts)
        h_tr1 = gROOT.FindObject('h_tr1')
        y_tr1 = h_tr1.GetMean()

        tree_data.Draw('tr2_cluster_y[0]>>h_tr2(64, 80, 195.2)', cuts)
        h_tr2 = gROOT.FindObject('h_tr2')
        y_tr2 = h_tr2.GetMean()

        tree_data.Draw('cal_cluster_y[0]>>h_cal(64, 80, 195.2)', cuts)
        h_cal = gROOT.FindObject('h_cal')
        y_cal = h_cal.GetMean()

        track = TGraphErrors(3, x, np.array([y_tr1, y_tr2, y_cal]))
        track.Fit('pol0', "Q")
        fit_func = track.GetFunction('pol0')
        print("Tr1 needs alignment for:", y_tr1 - fit_func.Eval(0))
        print("Tr2 needs alignment for:", y_tr2 - fit_func.Eval(0))
        print("Cal needs alignment for:", y_cal - fit_func.Eval(0))

    @classmethod
    def check_layer3_energy(cls):
        canvas = TCanvas('energy_2nd_cal_layer', 'energy_2nd_cal_layer', 1024, 768)
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
        h_mc.SetTitle('energy_2nd_cal_layer')
        gStyle.SetOptStat(1110)
        canvas.Write('energy_2nd_cal_layer')


class Calorimeter(Detector):

    def y(self):
        name = 'cal_y'
        canvas = TCanvas(name, name, 1024, 768)
        histos = []

        for i, tree in enumerate(self.trees):
            tree.Draw('cal_cluster_y>>h_y{}(64, 80, 195.2)'.format(i + 1))
            histos.append(gROOT.FindObject('h_y{}'.format(i + 1)))

        for i, tree in enumerate(self.trees):
            histos[i].Scale(1. / histos[i].GetEntries())
            histos[i].SetLineWidth(3)
            histos[i].SetLineColor(i + 1)
            histos[i].SetTitle('Beam Energy {} GeV'.format(i + 1))
            if i == 0:
                histos[i].GetXaxis().SetTitle('y, mm')
                histos[i].GetYaxis().SetTitle('Events, %')
                histos[i].Draw("histo")
            else:
                histos[i].Draw("histosame")

        canvas.BuildLegend()
        histos[0].SetTitle('Calorimeter')
        gStyle.SetOptStat(1110)
        canvas.Write()

    def energy(self):
        name = 'cal_energy'
        canvas = TCanvas(name, name, 1024, 768)
        histos = []

        for i, tree in enumerate(self.trees):
            tree.Draw('cal_cluster_energy>>h_energy{}(100, 0, 800)'.format(i + 1))
            histos.append(gROOT.FindObject('h_energy{}'.format(i + 1)))

        for i, tree in enumerate(self.trees):
            histos[i].Scale(1. / histos[i].GetEntries())
            histos[i].SetLineWidth(3)
            histos[i].SetLineColor(i + 1)
            histos[i].SetTitle('Beam energy {} GeV'.format(i + 1))
            if i == 0:
                histos[i].GetXaxis().SetTitle('E_{clst}, MIPs')
                histos[i].GetYaxis().SetTitle('Events, %')
                histos[i].Draw("histo")
            else:
                histos[i].Draw("histosame")

        canvas.BuildLegend()
        histos[0].SetTitle('Calorimeter')
        gStyle.SetOptStat(1110)
        canvas.Write()

    def energy_layer1(self):
        name = 'cal_energy_layer1'
        canvas = TCanvas(name, name, 1024, 768)
        histos = []

        for i, tree in enumerate(self.trees):
            tree.Draw('Sum$(cal_hit_energy*(cal_hit_layer == 2))>>h_energy{}(100, 0, 150)'.format(i + 1))
            histos.append(gROOT.FindObject('h_energy{}'.format(i + 1)))

        for i, tree in enumerate(self.trees):
            histos[i].Scale(1. / histos[i].GetEntries())
            histos[i].SetLineWidth(3)
            histos[i].SetLineColor(i + 1)
            histos[i].SetTitle('Beam energy {} GeV'.format(i + 1))
            if i == 0:
                histos[i].GetXaxis().SetTitle('E_{clst}, MIPs')
                histos[i].GetYaxis().SetTitle('Events, %')
                histos[i].Draw("histo")
            else:
                histos[i].Draw("histosame")

        canvas.BuildLegend()
        histos[0].SetTitle('Calorimeter layer 1 energy deposit')
        gStyle.SetOptStat(1110)
        canvas.Write()

    def n_hits_layer1(self):
        name = 'cal_n_hits_layer1'
        canvas = TCanvas(name, name, 1024, 768)
        histos = []

        for i, tree in enumerate(self.trees):
            tree.Draw('Sum$(cal_hit_layer == 2)>>h_n_hits{}(20, 0, 20)'.format(i + 1))
            histos.append(gROOT.FindObject('h_n_hits{}'.format(i + 1)))

        for i, tree in enumerate(self.trees):
            histos[i].Scale(1. / histos[i].GetEntries())
            histos[i].SetLineWidth(3)
            histos[i].SetLineColor(i + 1)
            histos[i].SetTitle('Beam energy {} GeV'.format(i + 1))
            if i == 0:
                histos[i].GetXaxis().SetTitle('N_{hits}')
                histos[i].GetYaxis().SetTitle('Events, %')
                histos[i].Draw("histo")
            else:
                histos[i].Draw("histosame")

        canvas.BuildLegend()
        histos[0].SetTitle('Calorimeter layer 1 n hits')
        gStyle.SetOptStat(1110)
        canvas.Write()

    def energy_1st_cal_layer_vs_beam_energy(self):
        name = 'energy_1st_cal_layer_vs_beam_energy'
        canvas = TCanvas(name, name, 1024, 768)
        canvas.cd()

        histo = TH2F(name, name, 5, 0.5, 5.5, 50, 0, 150)

        cuts = 'tr1_n_clusters == 1 && tr2_n_clusters == 1'
        self.tree_5gev.Draw('Sum$(cal_hit_energy*(cal_hit_layer == 2))>>h_5gev(50, 0, 150)', cuts)
        self.tree_4gev.Draw('Sum$(cal_hit_energy*(cal_hit_layer == 2))>>h_4gev(50, 0, 150)', cuts)
        self.tree_3gev.Draw('Sum$(cal_hit_energy*(cal_hit_layer == 2))>>h_3gev(50, 0, 150)', cuts)
        self.tree_2gev.Draw('Sum$(cal_hit_energy*(cal_hit_layer == 2))>>h_2gev(50, 0, 150)', cuts)
        self.tree_1gev.Draw('Sum$(cal_hit_energy*(cal_hit_layer == 2))>>h_1gev(50, 0, 150)', cuts)

        h_5gev = gROOT.FindObject('h_5gev')
        h_4gev = gROOT.FindObject('h_4gev')
        h_3gev = gROOT.FindObject('h_3gev')
        h_2gev = gROOT.FindObject('h_2gev')
        h_1gev = gROOT.FindObject('h_1gev')
        h_5gev.Scale(1. / self.tree_5gev.GetEntries())
        h_4gev.Scale(1. / self.tree_4gev.GetEntries())
        h_3gev.Scale(1. / self.tree_3gev.GetEntries())
        h_2gev.Scale(1. / self.tree_2gev.GetEntries())
        h_1gev.Scale(1. / self.tree_1gev.GetEntries())

        for i in range(1, 51):
            histo.Fill(5, h_5gev.GetBinCenter(i), h_5gev.GetBinContent(i))
            histo.Fill(4, h_4gev.GetBinCenter(i), h_4gev.GetBinContent(i))
            histo.Fill(3, h_3gev.GetBinCenter(i), h_3gev.GetBinContent(i))
            histo.Fill(2, h_2gev.GetBinCenter(i), h_2gev.GetBinContent(i))
            histo.Fill(1, h_1gev.GetBinCenter(i), h_1gev.GetBinContent(i))

        histo.Draw("COLZ")
        histo.GetXaxis().SetTitle("Beam energy, [GeV]")
        histo.GetYaxis().SetTitle("1st cal layer energy, [MIP]")
        histo.Write(name)

    def n_hits_1st_cal_layer_vs_beam_energy(self):
        name = 'n_hits_in_1st_cal_layer_vs_beam_energy'
        canvas = TCanvas(name, name, 1024, 768)
        canvas.cd()

        histo = TH2F(name, name, 5, 0.5, 5.5, 20, 0, 20)

        cuts = 'tr1_n_clusters == 1 && tr2_n_clusters == 1'
        self.tree_5gev.Draw('Sum$(cal_hit_layer == 2)>>h_5gev(20, 0, 20)', cuts)
        self.tree_4gev.Draw('Sum$(cal_hit_layer == 2)>>h_4gev(20, 0, 20)', cuts)
        self.tree_3gev.Draw('Sum$(cal_hit_layer == 2)>>h_3gev(20, 0, 20)', cuts)
        self.tree_2gev.Draw('Sum$(cal_hit_layer == 2)>>h_2gev(20, 0, 20)', cuts)
        self.tree_1gev.Draw('Sum$(cal_hit_layer == 2)>>h_1gev(20, 0, 20)', cuts)

        h_5gev = gROOT.FindObject('h_5gev')
        h_4gev = gROOT.FindObject('h_4gev')
        h_3gev = gROOT.FindObject('h_3gev')
        h_2gev = gROOT.FindObject('h_2gev')
        h_1gev = gROOT.FindObject('h_1gev')
        h_5gev.Scale(1. / self.tree_5gev.GetEntries())
        h_4gev.Scale(1. / self.tree_4gev.GetEntries())
        h_3gev.Scale(1. / self.tree_3gev.GetEntries())
        h_2gev.Scale(1. / self.tree_2gev.GetEntries())
        h_1gev.Scale(1. / self.tree_1gev.GetEntries())

        for i in range(1, 51):
            histo.Fill(5, h_5gev.GetBinCenter(i), h_5gev.GetBinContent(i))
            histo.Fill(4, h_4gev.GetBinCenter(i), h_4gev.GetBinContent(i))
            histo.Fill(3, h_3gev.GetBinCenter(i), h_3gev.GetBinContent(i))
            histo.Fill(2, h_2gev.GetBinCenter(i), h_2gev.GetBinContent(i))
            histo.Fill(1, h_1gev.GetBinCenter(i), h_1gev.GetBinContent(i))

        histo.Draw("COLZ")
        histo.GetXaxis().SetTitle("Beam energy, [GeV]")
        histo.GetYaxis().SetTitle("N_{hits} in 1st cal layer")
        histo.Write(name)


class Tracker(Detector):

    def __init__(self, tr_idx):
        self.tr_idx = tr_idx
        super().__init__()

    def n_clst_in_tr(self):
        '''Looks wheter 2 clusters in trackers appear more with increased energy'''
        name = 'n_clst_in_tr'
        canvas = TCanvas(name, name, 1024, 768)
        canvas.Divide(3, 1)

        energy_arr = [1., 2., 3., 4., 5.]

        events_ineff_arr = []
        events_ineff_err_arr = []
        for i, tree in enumerate(self.trees):
            tree.Draw(">>list_ineff{}".format(i), "tr{}_n_clusters <1".format(self.tr_idx))
            ev_list = gDirectory.Get("list_ineff{}".format(i))
            events_ineff_arr.append(ev_list.GetN() * 100 / tree.GetEntries())
            events_ineff_err_arr.append(ev_list.GetN()**0.5 * 100 / tree.GetEntries())

        gr_ineff = TGraphErrors(5, np.array(energy_arr), np.array(events_ineff_arr), np.zeros(5), np.array(events_ineff_err_arr))
        gr_ineff.SetTitle("Tr{} Ineficiency (0 hits)".format(self.tr_idx))

        canvas.cd(1)
        gr_ineff.Draw("APE")
        gr_ineff.GetXaxis().SetTitle("Beam energy, GeV")
        gr_ineff.GetYaxis().SetTitle("Events, %")

        events_arr = []
        events_err_arr = []
        for i, tree in enumerate(self.trees):
            tree.Draw(">>list{}".format(i), "tr{}_n_clusters == 1".format(self.tr_idx))
            ev_list = gDirectory.Get("list{}".format(i))
            events_arr.append(ev_list.GetN() * 100 / tree.GetEntries())
            events_err_arr.append(ev_list.GetN()**0.5 * 100 / tree.GetEntries())

        gr_ok = TGraphErrors(5, np.array(energy_arr), np.array(events_arr), np.zeros(5), np.array(events_err_arr))
        gr_ok.SetTitle("Tr{} is ok (1 hits)".format(self.tr_idx))

        canvas.cd(2)
        gr_ok.Draw("APE")
        gr_ok.GetXaxis().SetTitle("Beam energy, GeV")
        gr_ok.GetYaxis().SetTitle("Events, %")

        events_bs_arr = []
        events_bs_err_arr = []
        for i, tree in enumerate(self.trees):
            tree.Draw(">>list_bs{}".format(i), "tr{}_n_clusters > 1".format(self.tr_idx))
            ev_list = gDirectory.Get("list_bs{}".format(i))
            events_bs_arr.append(ev_list.GetN() * 100 / tree.GetEntries())
            events_bs_err_arr.append(ev_list.GetN()**0.5 * 100 / tree.GetEntries())

        gr_bs = TGraphErrors(5, np.array(energy_arr), np.array(events_bs_arr), np.zeros(5), np.array(events_bs_err_arr))
        gr_bs.SetTitle("Tr{} BackScattered (>1 hits)".format(self.tr_idx))

        canvas.cd(3)
        gr_bs.Draw("APE")
        gr_bs.GetXaxis().SetTitle("Beam energy, GeV")
        gr_bs.GetYaxis().SetTitle("Events, %")

        canvas.Write()

    def n_clst_in_tr_corr(self):
        '''Looks wheter 2 clusters in trackers appear more with increased energy'''
        name = 'tr{}_n_clst_corr_to_cal_n_hits'.format(self.tr_idx)
        canvas = TCanvas(name, name, 1024, 768)

        n_hits_cal_arr = [2.2913242, 3.1659892, 3.6743501, 4.0422381, 4.2624919]
        n_hits_cal_arr_err = []
        for point in n_hits_cal_arr:
            n_hits_cal_arr_err.append(point**0.5)

        events_bs_arr = []
        events_bs_err_arr = []
        for i, tree in enumerate(self.trees):
            tree.Draw(">>list_bs{}".format(i), "tr{}_n_clusters > 1".format(self.tr_idx))
            ev_list = gDirectory.Get("list_bs{}".format(i))
            events_bs_arr.append(ev_list.GetN() * 100 / tree.GetEntries())
            events_bs_err_arr.append(ev_list.GetN()**0.5 * 100 / tree.GetEntries())

        gr_bs = TGraphErrors(5, np.array(n_hits_cal_arr), np.array(events_bs_arr), np.array(n_hits_cal_arr_err), np.array(events_bs_err_arr))
        gr_bs.SetTitle("Correlation: N_clst_tr{} vs N_hits_layer1 ".format(self.tr_idx))

        gr_bs.Draw("APE")
        gr_bs.GetXaxis().SetTitle("Mean of N_{hits} in layer1 calorimeter")
        gr_bs.GetYaxis().SetTitle("Events with >1 cluster in tr{}, %".format(self.tr_idx))

        canvas.Write()

    def y(self):
        name = 'tr{}_y'.format(self.tr_idx)
        canvas = TCanvas(name, name, 1024, 768)
        histos = []

        for i, tree in enumerate(self.trees):
            tree.Draw('tr{}_cluster_y>>h_y{}(64, 80, 195.2)'.format(self.tr_idx, i + 1))
            histos.append(gROOT.FindObject('h_y{}'.format(i + 1)))

        for i, tree in enumerate(self.trees):
            histos[i].Scale(1. / histos[i].GetEntries())
            histos[i].SetLineWidth(3)
            histos[i].SetLineColor(i + 1)
            histos[i].SetTitle('Beam energy {} GeV'.format(i + 1))
            if i == 0:
                histos[i].GetXaxis().SetTitle('y, mm')
                histos[i].GetYaxis().SetTitle('Events, %')
                histos[i].Draw("histo")
            else:
                histos[i].Draw("histosame")

        canvas.BuildLegend()
        histos[0].SetTitle('Tracker{}'.format(self.tr_idx))
        gStyle.SetOptStat(1110)
        canvas.Write()

    def energy(self):
        name = 'tr{}_energy'.format(self.tr_idx)
        canvas = TCanvas(name, name, 1024, 768)
        histos = []

        for i, tree in enumerate(self.trees):
            tree.Draw('tr{}_cluster_energy>>h_energy{}(100, 0, 10)'.format(self.tr_idx, i + 1))
            histos.append(gROOT.FindObject('h_energy{}'.format(i + 1)))

        for i, tree in enumerate(self.trees):
            histos[i].Scale(1. / histos[i].GetEntries())
            histos[i].SetLineWidth(3)
            histos[i].SetLineColor(i + 1)
            histos[i].SetTitle('Beam energy {} GeV'.format(i + 1))
            if i == 0:
                histos[i].GetXaxis().SetTitle('E_{clst}, MIPs')
                histos[i].GetYaxis().SetTitle('Events, %')
                histos[i].Draw("histo")
            else:
                histos[i].Draw("histosame")

        canvas.BuildLegend()
        histos[0].SetTitle('Tracker{}'.format(self.tr_idx))
        gStyle.SetOptStat(1110)
        canvas.Write()

    def shower_distance_scan_energy(self):
        name = 'tr{}_shower_distance_scan_energy'.format(self.tr_idx)
        canvas = TCanvas(name, name, 1024, 768)
        canvas.cd()

        cuts = ''
        distance = 'tr{}_cluster_y - cal_cluster_y[0]'.format(self.tr_idx)
        self.tree_5gev.Draw(distance + '>>h_5gev(20, -70, 50)', cuts)
        self.tree_4gev.Draw(distance + '>>h_4gev(20, -70, 50)', cuts)
        self.tree_3gev.Draw(distance + '>>h_3gev(20, -70, 50)', cuts)
        self.tree_2gev.Draw(distance + '>>h_2gev(20, -70, 50)', cuts)
        self.tree_1gev.Draw(distance + '>>h_1gev(20, -70, 50)', cuts)

        h_5gev = gROOT.FindObject('h_5gev')
        h_4gev = gROOT.FindObject('h_4gev')
        h_3gev = gROOT.FindObject('h_3gev')
        h_2gev = gROOT.FindObject('h_2gev')
        h_1gev = gROOT.FindObject('h_1gev')
        h_5gev.Scale(1 / self.tree_5gev.GetEntries())
        h_4gev.Scale(1 / self.tree_4gev.GetEntries())
        h_3gev.Scale(1 / self.tree_3gev.GetEntries())
        h_2gev.Scale(1 / self.tree_2gev.GetEntries())
        h_1gev.Scale(1 / self.tree_1gev.GetEntries())

        histos = [h_1gev, h_2gev, h_3gev, h_4gev, h_5gev]

        for i, hist in enumerate(histos):
            if (i == 0):
                hist.Draw('histo')
                hist.SetTitle("Beam energy 1 GeV")
                hist.GetXaxis().SetTitle("d_{shower}, mm")
                hist.GetYaxis().SetTitle("Events, %")
                continue
            hist.Draw('histosame')
            hist.SetTitle("Beam energy {} GeV".format(i + 1))
            hist.SetLineColor(i + 1)

        canvas.SetLogy()
        canvas.BuildLegend()
        histos[0].SetTitle("Tracker{}".format(self.tr_idx))

        canvas.Write(name)


def main():

    cal = Calorimeter()

    cal.energy_1st_cal_layer_vs_beam_energy()
    cal.n_hits_1st_cal_layer_vs_beam_energy()
    input('Yaay I am finished :3')


gROOT.SetBatch(1)
gROOT.SetStyle('ATLAS')

main()
