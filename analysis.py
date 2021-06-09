import numpy as np
from ROOT import TFile, TCanvas, gROOT, THStack




def study_smearing():
    # Import data and MC files:
    f_data = TFile.Open("../trees_5gev_e/data_nn50.root", 'read')
    t_data = f_data.data
    t_data.AddFriend("data10=data", "../trees_5gev_e/data_nn10.root")
    t_data.AddFriend("data20=data", "../trees_5gev_e/data_nn20.root")
    t_data.AddFriend("data90=data", "../trees_5gev_e/data_nn90.root")
    t_data.AddFriend("data99=data", "../trees_5gev_e/data_nn99.root")

    f_mc = TFile.Open("../trees_5gev_e/lucas.root", 'read')
    t_mc = f_mc.lumical
    t_mc.AddFriend("lumical_noise=lumical", "../trees_5gev_e/lucas_geocuts_noise_test.root")

    c = TCanvas()
    hs = THStack("hs", "title")

    t_data.Draw("Sum$(data10.energy*(data10.layer == 0))>>h_data10(200, 0, 3)", "Sum$(data10.energy*(data10.layer == 0)) > 0")
    h_data10 = gROOT.FindObject("h_data10")
    h_data10.SetTitle("Data nn>0.1")
    h_data10.SetLineWidth(2)
    h_data10.SetLineColor(3)
    hs.Add(h_data10, "histo")

    t_data.Draw("Sum$(data20.energy*(data20.layer == 0))>>h_data20(200, 0, 3)", "Sum$(data20.energy*(data20.layer == 0)) > 0")
    h_data20 = gROOT.FindObject("h_data20")
    h_data20.SetTitle("Data nn>0.2")
    h_data20.SetLineWidth(2)
    h_data20.SetLineColor(4)
    hs.Add(h_data20, "histo")

    t_data.Draw("Sum$(energy*(layer == 0))>>h_data50(200, 0, 3)", "Sum$(energy*(layer == 0)) > 0")
    h_data50 = gROOT.FindObject("h_data50")
    h_data50.SetTitle("Data nn>0.5")
    h_data50.SetLineWidth(2)
    h_data50.SetLineColor(1)
    hs.Add(h_data50, "histo")

    t_data.Draw("Sum$(data90.energy*(data90.layer == 0))>>h_data90(200, 0, 3)", "Sum$(data90.energy*(data90.layer == 0)) > 0")
    h_data90 = gROOT.FindObject("h_data90")
    h_data90.SetTitle("Data nn>0.90")
    h_data90.SetLineWidth(2)
    h_data90.SetLineColor(5)
    hs.Add(h_data90, "histo")

    t_data.Draw("Sum$(data99.energy*(data99.layer == 0))>>h_data99(200, 0, 3)", "Sum$(data99.energy*(data99.layer == 0)) > 0")
    h_data99 = gROOT.FindObject("h_data99")
    h_data99.SetTitle("Data nn>0.99")
    h_data99.SetLineWidth(2)
    h_data99.SetLineColor(6)
    hs.Add(h_data99, "histo")

    # t_mc.Draw("Sum$(tr1_energy/0.0885)>>h_mc(200, 0, 3)", "Sum$(tr1_energy/0.0885) > 0")
    # h_mc = gROOT.FindObject("h_mc")
    # h_mc.SetTitle("MC no smearing")
    # h_mc.SetLineWidth(2)
    # h_mc.SetLineColor(2)
    # h_mc.Scale(t_data.GetEntries()/t_mc.GetEntries())
    # hs.Add(h_mc, "histo")

    # t_mc.Draw("Sum$(tr1_energy_smeared5/0.0885)>>h_mc5(200, 0, 3)", "Sum$(tr1_energy_smeared5/0.0885) > 0")
    # h_mc5 = gROOT.FindObject("h_mc5")
    # h_mc5.SetLineWidth(2)
    # h_mc5.SetLineColor(3)
    # h_mc5.SetMarkerStyle(0)
    # h_mc5.Scale(t_data.GetEntries()/t_mc.GetEntries())
    # hs.Add(h_mc5, "histo")

    # t_mc.Draw("Sum$(tr1_energy_smeared6/0.0885)>>h_mc6(200, 0, 3)", "Sum$(tr1_energy_smeared6/0.0885) > 0")
    # h_mc6 = gROOT.FindObject("h_mc6")
    # h_mc6.SetTitle("MC smearing #sigma=0.6*stdNoise")
    # h_mc6.SetLineWidth(2)
    # h_mc6.SetLineColor(2)
    # h_mc6.Scale(t_data.GetEntries()/t_mc.GetEntries())
    # hs.Add(h_mc6, "histo")

    t_mc.Draw("Sum$(tr1_energy_smeared7/0.0885)>>h_mc7(200, 0, 3)", "Sum$(tr1_energy_smeared7/0.0885) > 0")
    h_mc7 = gROOT.FindObject("h_mc7")
    h_mc7.SetTitle("MC smearing #sigma=0.7*stdNoise")
    h_mc7.SetLineWidth(2)
    h_mc7.SetLineColor(2)
    h_mc7.Scale(t_data.GetEntries()/t_mc.GetEntries())
    hs.Add(h_mc7, "histo")

    # t_mc.Draw("Sum$((tr1_energy_smeared10)/0.0885)>>h_mc10(200, 0, 3)", "Sum$((tr1_energy_smeared10)/0.0885) > 0")
    # h_mc10 = gROOT.FindObject("h_mc10")
    # h_mc10.SetTitle("MC smearing #sigma=stdNoise")
    # h_mc10.SetLineWidth(2)
    # h_mc10.SetLineColor(18)
    # h_mc10.Scale(t_data.GetEntries()/t_mc.GetEntries())
    # hs.Add(h_mc10, "histo")

    hs.Draw("nostack")
    hs.SetTitle("Deposited energy in the 1st tracker;E, [MIP]; N_{events}")

    c.SetGridx()
    c.SetGridy()
    c.BuildLegend()
    c.Update()
    input("wait")


def study_cal_eff():
    # Import data and MC files:
    f_data = TFile.Open("../trees_5gev_e/data_nn50.root", 'read')
    t_data = f_data.data

    f_mc = TFile.Open("../trees_5gev_e/lucas.root", 'read')
    t_mc = f_mc.lumical
    t_mc.AddFriend("lumical_noise=lumical", "../trees_5gev_e/lucas_geocuts_noise_cal.root")
    t_mc.AddFriend("lumical_eff=lumical", "../trees_5gev_e/lucas_geocuts_noise_cal_eff.root")

    c = TCanvas()
    hs = THStack("hs", "title")

    t_data.Draw("energy*(layer > 1)>>h_data50(200, 0, 50)", "energy*(layer > 1) > 0")
    h_data50 = gROOT.FindObject("h_data50")
    h_data50.SetTitle("Data")
    h_data50.SetLineWidth(2)
    h_data50.SetLineColor(1)
    hs.Add(h_data50, "histo")

    t_mc.Draw("cal_energy_new/0.0885>>h_mc(200, 0, 50)", "cal_energy_new/0.0885 > 0","", 100000)
    h_mc = gROOT.FindObject("h_mc")
    h_mc.SetTitle("MC smearing and cal eff")
    h_mc.SetLineWidth(2)
    h_mc.SetLineColor(3)
    h_mc.Scale(t_data.GetEntries()/100000)#t_mc.GetEntries())
    hs.Add(h_mc, "histo")

    t_mc.Draw("cal_energy_smeared7/0.0885>>h_mc7(200, 0, 50)", "cal_energy_smeared7/0.0885 > 0", "", 100000)
    h_mc7 = gROOT.FindObject("h_mc7")
    h_mc7.SetTitle("MC smearing=0.7*#sigma")
    h_mc7.SetLineWidth(2)
    h_mc7.SetLineColor(2)
    h_mc7.Scale(t_data.GetEntries()/100000)#t_mc.GetEntries())
    hs.Add(h_mc7, "histo")
    print(h_mc7.GetEntries())
    hs.Draw("nostack")
    hs.SetTitle("Deposited energy in the pads of calorimeter;E, [MIP]; N_{hits}")

    c.SetGridx()
    c.SetGridy()
    c.BuildLegend()
    c.Update()
    input("wait")


def study_trigger():
    # Import data and MC files:
    f_data = TFile.Open("../trees_5gev_e/.root", 'read')
    t_data = f_data.data

    f_mc = TFile.Open("../trees_5gev_e/lucas.root", 'read')
    t_mc = f_mc.lumical
    t_mc.AddFriend("lumical_tr_clusters=lumical", "../trees_5gev_e/.root")

    c = TCanvas()

    t_data.Draw("tr1_n_clusters>>h_data(6, 0, 6)")
    h_data = gROOT.FindObject("h_data")
    h_data.SetTitle("Data")
    h_data.SetLineWidth(2)
    h_data.SetLineColor(1)

    t_mc.Draw("tr1_n_hits>>h_mc(6, 0, 6)")
    h_mc = gROOT.FindObject("h_mc")
    h_mc.SetTitle("MC")
    h_mc.SetLineWidth(2)
    h_mc.SetLineColor(2)
    h_mc.Scale(t_data.GetEntries()/t_mc.GetEntries())

    t_mc.Draw("tr1_n_clusters>>h_mc2(6, 0, 6)", "n_triggers == 3")
    h_mc2 = gROOT.FindObject("h_mc2")
    h_mc2.SetTitle("MC w/ trigger")
    h_mc2.SetLineWidth(2)
    h_mc2.SetLineColor(4)
    h_mc2.Scale(t_data.GetEntries()/t_mc.GetEntries())


    hs = THStack("hs", "title")
    hs.Add(h_data, "histo")
    hs.Add(h_mc, "histo")
    hs.Add(h_mc2, "histo")
    hs.Draw("nostack")
    hs.SetTitle("Number of clusters in the 1st tracker;N_{particles}; N_{events}")

    c.SetGridx()
    c.SetGridy()
    c.BuildLegend()
    c.Update()
    input("wait")


def study_geo_cuts():
    # Import data and MC files:
    f_data = TFile.Open("../trees_5gev_e/data_nn50.root", 'read')
    t_data = f_data.data

    f_mc = TFile.Open("../trees_5gev_e/lucas.root", 'read')
    t_mc = f_mc.lumical
    t_mc.AddFriend("lumical_tr_clusters=lumical", "../trees/lucas_5gev_e_tr_clusters.root")

    c = TCanvas("tr1_mc_hit_map")

    t_data.Draw("cal_hit_pad:cal_hit_sector>>h_data(4, 0, 4, 64, 0, 64)", "cal_hit_layer == 6", "colz")
    h_data = gROOT.FindObject("h_data")
    h_data.SetTitle("Data: hit map in tracker 2; sector, [0-3]; pad, [0-63]")
    h_data.SetStats(0)

    # t_mc.Draw("tr1_pad:tr1_sector>>h_mc(4, 0, 4, 64, 0, 64)", "", "colz")
    # h_mc = gROOT.FindObject("h_mc")
    # h_mc.SetTitle("MC: hit map in tracker 1; sector, [0-3]; pad, [0-63]")
    # h_mc.SetStats(0)

    c.SetLogz()
    c.Update()
    input("wait")


# study_smearing()
study_cal_eff()
# study_trigger()
