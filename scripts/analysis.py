from ROOT import TFile, gROOT, TGraphErrors, TH1F, TGraph, TCanvas, TPad, TF1, gStyle
import numpy as np


def main():
    file_mc = TFile.Open("../extracted_root_files/extracted_mc_5GeV.root")
    tree_mc = file_mc.mc

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

    tree_mc.Draw('tr1_hit_y - cal_cluster_y[0]>>h_mc(64, -32., 32)')
    histo = gROOT.FindObject('h_mc')
    histo.Scale(1. / tree_mc.GetEntries())
    histo.SetLineWidth(3)
    histo.SetLineColor(6)
    histo.SetFillColor(6)
    histo.SetTitle('Monte Carlo 5 GeV')

    for i, tree in enumerate(trees):
        if i == 0:
            tree.Draw('tr1_hit_y-cal_cluster_y[0]>>h_energy{}(64, -32., 32)'.format(i + 1), "", "histosame")
        else:
            tree.Draw('tr1_hit_y-cal_cluster_y[0]>>h_energy{}(64, -32., 32)'.format(i + 1), "", "histosame")

        histo = gROOT.FindObject('h_energy{}'.format(i + 1))
        histo.Scale(1. / tree.GetEntries())
        histo.SetLineWidth(3)
        histo.SetLineColor(i + 1)
        histo.SetTitle('Beam energy {} GeV'.format(i + 1))

    canvas.BuildLegend()
    canvas.Write()

    # tree_mc.Draw("hit_pad>>h_all(64, 0, 64)", "hit_energy*(hit_layer == 1)", "HISTO")
    # h_all = gROOT.FindObject("h_all")
    # h_all.SetLineWidth(3)
    # h_all.SetTitle("All hits")
    # h_all.Write()

    # tree_mc.Draw("hit_pad>>h_bs(64, 0, 64)", "hit_energy*(hit_layer == 1 && hit_bs == 1)", "samehisto")
    # h_bs = gROOT.FindObject("h_bs")
    # h_bs.SetLineColor(2)
    # h_bs.SetLineWidth(3)
    # h_bs.SetTitle("Back-scattered hits")
    # h_bs.Write()

    # tree_mc.Draw("hit_pad>>h_dir(64, 0, 64)", "hit_energy*(hit_layer == 1 && hit_bs == 0)", "samehisto")
    # h_dir = gROOT.FindObject("h_dir")
    # h_dir.SetLineColor(4)
    # h_dir.SetLineWidth(3)
    # h_dir.SetTitle("Direct hits")
    # h_dir.Write()

    input('Yaay I am finished :3')


gROOT.SetBatch(0)
gROOT.SetStyle('ATLAS')

main()
