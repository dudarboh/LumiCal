from ROOT import TFile, gROOT, gDirectory, TGraphErrors, TCanvas, gStyle, nullptr, TH2F
import numpy as np


class Detector:
    file_5gev = TFile.Open('./extracted_trees/without_1.4_energy_cut/extracted_5_gev_energy_scan_1.root', 'read')
    tree_5gev = file_5gev.data

    file_4gev = TFile.Open('./extracted_trees/without_1.4_energy_cut/extracted_4_gev_energy_scan_1.root', 'read')
    tree_4gev = file_4gev.data

    file_3gev = TFile.Open('./extracted_trees/without_1.4_energy_cut/extracted_3_gev_energy_scan_1.root', 'read')
    tree_3gev = file_3gev.data

    file_2gev = TFile.Open('./extracted_trees/without_1.4_energy_cut/extracted_2_gev_energy_scan_1.root', 'read')
    tree_2gev = file_2gev.data

    file_1gev = TFile.Open('./extracted_trees/without_1.4_energy_cut/extracted_1_gev_energy_scan_1.root', 'read')
    tree_1gev = file_1gev.data

    trees = [tree_1gev, tree_2gev, tree_3gev, tree_4gev, tree_5gev]

    output_file = TFile('results.root', 'recreate')

    def __init__(self):
        self.output_file.cd()
        for i, tree in enumerate(self.trees):
            print(i + 1, "GeV,", "N_ev:", tree.GetEntries())


class Calorimeter(Detector):

    def variable_vs_layer(self):
        '''Change expression in TTree.Draw() to plot what you want
            None of the parameters couldnt see back-scattering in 1st (3rd) calorimeter layer
        '''
        name = 'std_energy_vs_layer'

        layers = np.array([2., 3., 4., 5., 6.])
        var = np.zeros((5, 5))
        error_var = np.zeros((5, 5))

        graphs = []

        for i, tree in enumerate(self.trees):
            for j, l in enumerate(layers):
                tree.Draw("cal_hit_y-cal_cluster_y[0]>>h_temp", "cal_hit_energy*(cal_hit_layer == {})".format(l))
                # tree.Draw("cal_hit_y-cal_cluster_y[0]>>h_temp", "(cal_hit_layer == {})".format(l))
                h_temp = gROOT.FindObject("h_temp")
                if i == 2:
                    h_temp.Scale(100.)
                var[i, j] = h_temp.GetStdDev()
                error_var[i, j] = h_temp.GetStdDevError()

            graphs.append(TGraphErrors(5, layers, var[i, :], nullptr, error_var[i, :]))
            graphs[i].SetTitle("%d GeV beam" % (i + 1))
            graphs[i].GetXaxis().SetTitle("layers")
            graphs[i].GetYaxis().SetTitle("StdDev, [mm]")
            graphs[i].SetMarkerStyle(20 + i)
            graphs[i].SetMarkerColor(i + 1)
            graphs[i].SetLineColor(i + 1)

        canvas = TCanvas(name, name, 1024, 768)
        canvas.cd()
        for i, graph in enumerate(graphs):
            if i == 0:
                graph.Draw('APL')
            else:
                graph.Draw('PLsame')

        canvas.BuildLegend()
        graphs[0].SetTitle("StdDev energy vs layer")
        gStyle.SetOptStat(1110)
        canvas.Write()

    # 5 histos on the pad
    def histo_vs_beam_energy(self):
        name = 'cal_energy_distr_layer1'
        canvas = TCanvas(name, name, 1024, 768)
        histos = []

        for i, tree in enumerate(self.trees):
            tree.Draw('cal_hit_y-cal_cluster_y[0]>>h_energy{}(6000, -100, 100)'.format(i + 1), "cal_hit_energy*(cal_hit_layer == 2)")
            histos.append(gROOT.FindObject('h_energy{}'.format(i + 1)))

        for i, histo in enumerate(histos):
            histo.Scale(1. / histo.GetEntries())
            histo.SetLineWidth(3)
            histo.SetLineColor(i + 1)
            histo.SetTitle('Beam energy {} GeV'.format(i + 1))
            if i == 0:
                histo.GetXaxis().SetTitle('y_{hit}-cluster_{main}, [mm]')
                histo.GetYaxis().SetTitle('Events, %')
                histo.SetMaximum(1.)
                histo.SetMinimum(0.)
                histo.Draw("histo")
            else:
                histo.Draw("histosame")

        canvas.BuildLegend()
        histos[0].SetTitle('Calorimeter layer 1 energy spectrum')
        gStyle.SetOptStat(1110)
        canvas.Write()


class Tracker(Detector):

    def __init__(self, tr_idx):
        self.tr_idx = tr_idx
        super().__init__()

    def histo_vs_beam_energy(self):
        name = 'tr{}_energy_spectrum'.format(self.tr_idx)
        canvas = TCanvas(name, name, 1024, 768)
        histos = []

        for i, tree in enumerate(self.trees):
            tree.Draw('tr{}_hit_energy>>h_energy{}(200, 0, 2.5)'.format(self.tr_idx, i + 1))
            histos.append(gROOT.FindObject('h_energy{}'.format(i + 1)))

        for i, tree in enumerate(self.trees):
            histos[i].Scale(1. / tree.GetEntries())
            histos[i].SetLineWidth(3)
            histos[i].SetLineColor(i + 1)
            histos[i].SetTitle('Beam energy {} GeV'.format(i + 1))
            if i == 0:
                histos[i].GetXaxis().SetTitle('E_{hit}, MIPs')
                histos[i].GetYaxis().SetTitle('Events, %')
                histos[i].SetMaximum(1.)
                histos[i].SetMinimum(0.)
                histos[i].Draw("histo")
            else:
                histos[i].Draw("histosame")

        canvas.BuildLegend()
        histos[0].SetTitle('Tracker{} layer 1 energy spectrum'.format(self.tr_idx))
        gStyle.SetOptStat(1110)
        canvas.Write()

    def shower_distance(self):
        name = 'tr{}_shower_distance_vs_beam_energy'.format(self.tr_idx)
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
    cal.variable_vs_layer()
    cal.histo_vs_beam_energy()

    # tr2 = Tracker(2)
    # tr2.histo_vs_beam_energy()

    # tr1 = Tracker(1)
    # tr1.histo_vs_beam_energy()

    input('Yaay I am finished :3')


gROOT.SetBatch(1)
gROOT.SetStyle('ATLAS')

main()
