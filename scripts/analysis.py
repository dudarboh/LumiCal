from ROOT import TFile, gROOT, TGraphErrors, TH1F, TGraph, TCanvas, TPad, TF1, gStyle
import numpy as np


def main():
    file_mc = TFile.Open("../extracted_root_files/extracted_mc_5GeV.root", 'read')
    tree_mc = file_mc.mc

    file_raw_mc = TFile.Open("../mc_root_files/LUCAS_TB16_5GeV.root", 'read')
    tree_raw_mc = file_raw_mc.LumiCal

    # file_mc_itamar = TFile.Open("/home/FoxWise/Documents/FCAL/analysis/mc_root_files/Itamar_mc_old/T16NST5G_29_12-10_35outputfile.root")
    # tree_mc_itamar = file_mc_itamar.Lcal

    file_5gev = TFile.Open('../extracted_root_files/with_1.4_energy_cut/extracted_5_gev_energy_scan_1.root', 'read')
    tree_5gev = file_5gev.data

    file_4gev = TFile.Open('../extracted_root_files/with_1.4_energy_cut/extracted_4_gev_energy_scan_1.root', 'read')
    tree_4gev = file_4gev.data

    file_3gev = TFile.Open('../extracted_root_files/with_1.4_energy_cut/extracted_3_gev_energy_scan_1.root', 'read')
    tree_3gev = file_3gev.data

    file_2gev = TFile.Open('../extracted_root_files/with_1.4_energy_cut/extracted_2_gev_energy_scan_1.root', 'read')
    tree_2gev = file_2gev.data

    file_1gev = TFile.Open('../extracted_root_files/with_1.4_energy_cut/extracted_1_gev_energy_scan_1.root', 'read')
    tree_1gev = file_1gev.data

    trees = [tree_1gev, tree_2gev, tree_3gev, tree_4gev, tree_5gev]
    # trees = [tree_5gev]
    #trees = [tree_1gev, tree_3gev, tree_5gev]

    output_file = TFile.Open("results.root", "RECREATE")
    output_file.cd()

    canvas = TCanvas("canvas", "title", 1024, 768)
    canvas.cd()

    # tree_mc.Draw('tr1_hit_y - cal_cluster_y[0]>>h_mc(64, -32., 32)')
    # histo = gROOT.FindObject('h_mc')
    # histo.Scale(1. / tree_mc.GetEntries())
    # histo.SetLineWidth(3)
    # histo.SetLineColor(6)
    # histo.SetFillColor(6)
    # histo.SetTitle('Monte Carlo 5 GeV')

    # for i, tree in enumerate(trees):
    #     if i == 0:
    #         tree.Draw('tr1_hit_y-cal_cluster_y[0]>>h_energy{}(64, -32., 32)'.format(i + 1), "", "histosame")
    #     else:
    #         tree.Draw('tr1_hit_y-cal_cluster_y[0]>>h_energy{}(64, -32., 32)'.format(i + 1), "", "histosame")

    #     histo = gROOT.FindObject('h_energy{}'.format(i + 1))
    #     histo.Scale(1. / tree.GetEntries())
    #     histo.SetLineWidth(3)
    #     histo.SetLineColor(i + 1)
    #     histo.SetTitle('Beam energy {} GeV'.format(i + 1))

    # I keep it here to mention after whether I did mistake here or not
    # h1 = TH1F('h1', "All hits", 100, -60., 40.)
    # h2 = TH1F('h2', "Back-scattered hits", 100, -60., 40.)
    # h3 = TH1F('h3', "Direct hits", 100, -60., 40.)
    # for idx, event in enumerate(tree_raw_mc):
    #     if idx % 1000 == 0:
    #         print(idx)

    #     hit_n = event.hit_n
    #     hit_layer = event.hit_layer
    #     hit_pad = event.hit_pad
    #     hit_bs = event.hit_bs
    #     hit_energy = event.hit_energy
    #     tot_energy = sum([hit_energy[i]/0.0885 if hit_layer[i] > 20 else 0. for i in range(hit_n)])
    #     weights = [max(0, 3.4 + np.log(hit_energy[i]/0.0885 / tot_energy)) if hit_layer[i] > 20 else 0. for i in range(hit_n)]
    #     pos = sum([weights[i] * (hit_pad[i] * 1.8 + 80. + 0.9) if hit_layer[i] > 20 else 0. for i in range(hit_n)])
    #     skipped_events = 0

    #     if sum(weights) != 0:
    #         pos /= sum(weights)
    #     else:
    #         skipped_events += 1
    #         print("I skipped:", skipped_events, "Events")

    #     for i in range(hit_n):
    #         if hit_layer[i] == 6:
    #             h1.Fill((hit_pad[i] * 1.8 + 80. + 0.9) - pos, hit_energy[i]/0.0885)

    #     for i in range(hit_n):
    #         if hit_layer[i] == 6 and hit_bs[i] == 1:
    #             h2.Fill((hit_pad[i] * 1.8 + 80. + 0.9) - pos, hit_energy[i]/0.0885)

    #     for i in range(hit_n):
    #         if hit_layer[i] == 6 and hit_bs[i] == 0:
    #             h3.Fill((hit_pad[i] * 1.8 + 80. + 0.9) - pos, hit_energy[i]/0.0885)

    # h1.Draw("histo")
    # h1.SetTitle("All hits")
    # h1.GetXaxis().SetTitle("y_{hit} - y_{main cluster}, mm")
    # h1.GetYaxis().SetTitle("Number of Hits * E_{hit}")
    # h1.SetLineColor(1)
    # h1.SetLineWidth(3)
    # h2.Draw("histosame")
    # h2.SetLineColor(2)
    # h2.SetLineWidth(3)
    # h3.Draw("histosame")
    # h3.SetLineColor(4)
    # h3.SetLineWidth(3)

    tree_mc.Draw("tr2_hit_sector>>h_all(4, 1, 5)", "", "HISTO")
    h_all = gROOT.FindObject("h_all")
    h_all.SetLineColor(1)
    h_all.SetLineWidth(3)
    h_all.SetTitle("All hits")
    h_all.GetYaxis().SetTitle("Number of Hits")
    h_all.GetXaxis().SetTitle("y_{hit} - y_{main clst}, mm")
    h_all.Write()

    tree_mc.Draw("tr2_hit_sector>>h_bs(4, 1, 5)", "(tr2_hit_bs == 1)", "samehisto")
    h_bs = gROOT.FindObject("h_bs")
    h_bs.SetLineColor(2)
    h_bs.SetLineWidth(3)
    h_bs.SetTitle("Back-scattered hits")
    h_bs.Write()

    tree_mc.Draw("tr2_hit_sector>>h_dir(4, 1, 5)", "(tr2_hit_bs == 0)", "samehisto")
    h_dir = gROOT.FindObject("h_dir")
    h_dir.SetLineColor(4)
    h_dir.SetLineWidth(3)
    h_dir.SetTitle("Direct hits")
    h_dir.Write()

    canvas.BuildLegend()
    canvas.Write()

    input('Yaay I am finished :3')


gROOT.SetBatch(0)
gROOT.SetStyle('ATLAS')

main()
