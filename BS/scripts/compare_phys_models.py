from ROOT import TFile, gROOT, TGraphErrors, TH1F, TGraph, TCanvas, TPad, TF1, gStyle, TF1, nullptr
import numpy as np

gROOT.SetBatch(1)
gROOT.SetStyle('ATLAS')

# Data trees
f_data_5gev = TFile.Open("../extracted_root_files/extracted_data_5gev.root", 'read')
f_data_4gev = TFile.Open('../extracted_root_files/extracted_data_4gev.root', 'read')
f_data_3gev = TFile.Open('../extracted_root_files/extracted_data_3gev.root', 'read')
f_data_2gev = TFile.Open('../extracted_root_files/extracted_data_2gev.root', 'read')
f_data_1gev = TFile.Open('../extracted_root_files/extracted_data_1gev.root', 'read')
t_data_5gev = f_data_5gev.data
t_data_4gev = f_data_4gev.data
t_data_3gev = f_data_3gev.data
t_data_2gev = f_data_2gev.data
t_data_1gev = f_data_1gev.data

# MC default physics model trees
f_mc_def_5gev = TFile.Open("../extracted_root_files/test_phys_models/extracted_lucas_5gev_def_phys.root", 'read')
f_mc_def_4gev = TFile.Open("../extracted_root_files/test_phys_models/extracted_lucas_4gev_def_phys.root", 'read')
f_mc_def_3gev = TFile.Open("../extracted_root_files/test_phys_models/extracted_lucas_3gev_def_phys.root", 'read')
f_mc_def_2gev = TFile.Open("../extracted_root_files/test_phys_models/extracted_lucas_2gev_def_phys.root", 'read')
f_mc_def_1gev = TFile.Open("../extracted_root_files/test_phys_models/extracted_lucas_1gev_def_phys.root", 'read')
t_mc_def_5gev = f_mc_def_5gev.mc
t_mc_def_4gev = f_mc_def_4gev.mc
t_mc_def_3gev = f_mc_def_3gev.mc
t_mc_def_2gev = f_mc_def_2gev.mc
t_mc_def_1gev = f_mc_def_1gev.mc

# MC Option3 EM physics model trees
f_mc_emy_5gev = TFile.Open("../extracted_root_files/test_phys_models/extracted_lucas_5gev_emy_phys.root", 'read')
f_mc_emy_4gev = TFile.Open("../extracted_root_files/test_phys_models/extracted_lucas_4gev_emy_phys.root", 'read')
f_mc_emy_3gev = TFile.Open("../extracted_root_files/test_phys_models/extracted_lucas_3gev_emy_phys.root", 'read')
f_mc_emy_2gev = TFile.Open("../extracted_root_files/test_phys_models/extracted_lucas_2gev_emy_phys.root", 'read')
f_mc_emy_1gev = TFile.Open("../extracted_root_files/test_phys_models/extracted_lucas_1gev_emy_phys.root", 'read')
t_mc_emy_5gev = f_mc_emy_5gev.mc
t_mc_emy_4gev = f_mc_emy_4gev.mc
t_mc_emy_3gev = f_mc_emy_3gev.mc
t_mc_emy_2gev = f_mc_emy_2gev.mc
t_mc_emy_1gev = f_mc_emy_1gev.mc

# MC Option4 EM physics model trees
f_mc_emz_5gev = TFile.Open("../extracted_root_files/test_phys_models/extracted_lucas_5gev_emz_phys.root", 'read')
f_mc_emz_4gev = TFile.Open("../extracted_root_files/test_phys_models/extracted_lucas_4gev_emz_phys.root", 'read')
f_mc_emz_3gev = TFile.Open("../extracted_root_files/test_phys_models/extracted_lucas_3gev_emz_phys.root", 'read')
f_mc_emz_2gev = TFile.Open("../extracted_root_files/test_phys_models/extracted_lucas_2gev_emz_phys.root", 'read')
f_mc_emz_1gev = TFile.Open("../extracted_root_files/test_phys_models/extracted_lucas_1gev_emz_phys.root", 'read')
t_mc_emz_5gev = f_mc_emz_5gev.mc
t_mc_emz_4gev = f_mc_emz_4gev.mc
t_mc_emz_3gev = f_mc_emz_3gev.mc
t_mc_emz_2gev = f_mc_emz_2gev.mc
t_mc_emz_1gev = f_mc_emz_1gev.mc

output_file = TFile.Open("./results.root", "RECREATE")
output_file.cd()

list_1gev = [t_data_1gev, t_mc_def_1gev, t_mc_emy_1gev, t_mc_emz_1gev]
list_2gev = [t_data_2gev, t_mc_def_2gev, t_mc_emy_2gev, t_mc_emz_2gev]
list_3gev = [t_data_3gev, t_mc_def_3gev, t_mc_emy_3gev, t_mc_emz_3gev]
list_4gev = [t_data_4gev, t_mc_def_4gev, t_mc_emy_4gev, t_mc_emz_4gev]
list_5gev = [t_data_5gev, t_mc_def_5gev, t_mc_emy_5gev, t_mc_emz_5gev]

list_energies = [list_1gev, list_2gev, list_3gev, list_4gev, list_5gev]

tr2_energy_podgonian = 1.08


def check_alignment(tree):
    x = np.array([0 * 4.5 + 2.25, 5 * 4.5 + 2.25, 23 * 4.5 + 2.25])

    cuts = "tr1_n_clusters == 1 && tr2_n_clusters == 1 && cal_n_clusters == 1"

    tree.Draw('tr1_cluster_y[0]>>h_tr1(64, 80, 195.2)', cuts)
    h_tr1 = gROOT.FindObject('h_tr1')
    y_tr1 = h_tr1.GetMean()

    tree.Draw('tr2_cluster_y[0]>>h_tr2(64, 80, 195.2)', cuts)
    h_tr2 = gROOT.FindObject('h_tr2')
    y_tr2 = h_tr2.GetMean()

    tree.Draw('cal_cluster_y[0]>>h_cal(64, 80, 195.2)', cuts)
    h_cal = gROOT.FindObject('h_cal')
    y_cal = h_cal.GetMean()

    track = TGraphErrors(3, x, np.array([y_tr1, y_tr2, y_cal]))
    track.Fit('pol0', "Q")
    fit_func = track.GetFunction('pol0')
    print("Tr1 needs alignment for:", y_tr1 - fit_func.Eval(0))
    print("Tr2 needs alignment for:", y_tr2 - fit_func.Eval(0))
    print("Cal needs alignment for:", y_cal - fit_func.Eval(0))


def control_plot(trees, variable="", binning="", selection="", canvas_name="name"):
    canvas = TCanvas(canvas_name, canvas_name, 1024, 768)

    titles = ["data", "Default MC", "Opt3 MC", "Opt4 MC"]
    colors = [1, 2, 4, 5, 6]
    histo_names = ["h_data", "h_mc_def", "h_mc_emy", "h_mc_emz"]
    for i, tree in enumerate(trees):
        if i == 0:
            if variable == "tr2_hit_energy":
                tree.Draw(variable + "*1.08" + ">>" + histo_names[i] + binning, selection, "histo")
            else:
                tree.Draw(variable + ">>" + histo_names[i] + binning, selection, "histo")

            histo = gROOT.FindObject(histo_names[i])
            histo.GetYaxis().SetTitle("#frac{N_{hits}}{N_{events}")
            histo.GetXaxis().SetTitle(variable)
        else:
            tree.Draw(variable + ">>" + histo_names[i] + binning, selection, "histosame")

        histo = gROOT.FindObject(histo_names[i])
        histo.Scale(1. / tree.GetEntries())
        histo.SetLineColor(colors[i])
        histo.SetTitle(titles[i])
        histo.SetLineWidth(3)

    canvas.BuildLegend()

    canvas.Write(canvas_name)


def avr_tr_energy(tr_idx):
    canvas = TCanvas("tr{}_energy_avr".format(tr_idx), "title", 1024, 768)
    x = [1., 2., 3., 4., 5.]
    y_data = []
    y_mc_def = []
    y_mc_emy = []
    y_mc_emz = []

    for tree in list_energies:
        tree[0].Draw('Sum$(tr{}_hit_energy)>>h_temp'.format(tr_idx))
        h_temp = gROOT.FindObject('h_temp')
        y_data.append(h_temp.GetMean())

        tree[1].Draw('Sum$(tr{}_hit_energy)>>h_temp'.format(tr_idx))
        h_temp = gROOT.FindObject('h_temp')
        y_mc_def.append(h_temp.GetMean())

        tree[2].Draw('Sum$(tr{}_hit_energy)>>h_temp'.format(tr_idx))
        h_temp = gROOT.FindObject('h_temp')
        y_mc_emy.append(h_temp.GetMean())

        tree[3].Draw('Sum$(tr{}_hit_energy)>>h_temp'.format(tr_idx))
        h_temp = gROOT.FindObject('h_temp')
        y_mc_emz.append(h_temp.GetMean())

    gr_data = TGraphErrors(5, np.array(x), np.array(y_data), nullptr, nullptr)
    gr_data.SetTitle("Data")
    gr_data.SetMarkerStyle(20)
    gr_data.SetMarkerColor(1)
    gr_data.Draw("AP")

    gr_mc_def = TGraphErrors(5, np.array(x), np.array(y_mc_def), nullptr, nullptr)
    gr_mc_def.SetTitle("Default MC")
    gr_mc_def.SetMarkerStyle(1)
    gr_mc_def.SetLineColor(2)
    gr_mc_def.Draw("Lsame")

    gr_mc_emy = TGraphErrors(5, np.array(x), np.array(y_mc_emy), nullptr, nullptr)
    gr_mc_emy.SetTitle("Opt3 MC")
    gr_mc_emy.SetMarkerStyle(1)
    gr_mc_emy.SetLineColor(4)
    gr_mc_emy.Draw("Lsame")

    gr_mc_emz = TGraphErrors(5, np.array(x), np.array(y_mc_emz), nullptr, nullptr)
    gr_mc_emz.SetTitle("Opt4 MC")
    gr_mc_emz.SetMarkerStyle(1)
    gr_mc_emz.SetLineColor(5)
    gr_mc_emz.Draw("Lsame")

    canvas.BuildLegend()
    canvas.Write("tr{}_energy_avr".format(tr_idx))
    input("stop")


# # Slide number 1
# control_plot(list_5gev, "tr1_hit_energy", "(200, 0, 5)", "", "tr1_energy_5gev")
# control_plot(list_5gev, "tr2_hit_energy", "(200, 0, 5)", "", "tr2_energy_5gev")
# control_plot(list_5gev, "tr1_hit_y-cal_cluster_y[0]", "(200, -60, 45)", "", "tr1_y_dist_5gev")
# control_plot(list_5gev, "tr2_hit_y-cal_cluster_y[0]", "(200, -60, 45)", "", "tr2_y_dist_5gev")
# # Slide number 2
# control_plot(list_5gev, "cal_hit_energy", "(200, 0, 100)", "", "cal_energy_5gev")

# # Slide number 3
# control_plot(list_4gev, "tr1_hit_energy", "(200, 0, 5)", "", "tr1_energy_4gev")
# control_plot(list_4gev, "tr2_hit_energy", "(200, 0, 5)", "", "tr2_energy_4gev")
# control_plot(list_4gev, "tr1_hit_y-cal_cluster_y[0]", "(200, -60, 45)", "", "tr1_y_dist_4gev")
# control_plot(list_4gev, "tr2_hit_y-cal_cluster_y[0]", "(200, -60, 45)", "", "tr2_y_dist_4gev")
# # Slide number 4
# control_plot(list_4gev, "cal_hit_energy", "(200, 0, 100)", "", "cal_energy_4gev")

# # Slide number 5
# control_plot(list_3gev, "tr1_hit_energy", "(200, 0, 5)", "", "tr1_energy_3gev")
# control_plot(list_3gev, "tr2_hit_energy", "(200, 0, 5)", "", "tr2_energy_3gev")
# control_plot(list_3gev, "tr1_hit_y-cal_cluster_y[0]", "(200, -60, 45)", "", "tr1_y_dist_3gev")
# control_plot(list_3gev, "tr2_hit_y-cal_cluster_y[0]", "(200, -60, 45)", "", "tr2_y_dist_3gev")
# # Slide number 6
# control_plot(list_3gev, "cal_hit_energy", "(200, 0, 100)", "", "cal_energy_3gev")

# # Slide number 7
# control_plot(list_2gev, "tr1_hit_energy", "(200, 0, 5)", "", "tr1_energy_2gev")
# control_plot(list_2gev, "tr2_hit_energy", "(200, 0, 5)", "", "tr2_energy_2gev")
# control_plot(list_2gev, "tr1_hit_y-cal_cluster_y[0]", "(200, -60, 45)", "", "tr1_y_dist_2gev")
# control_plot(list_2gev, "tr2_hit_y-cal_cluster_y[0]", "(200, -60, 45)", "", "tr2_y_dist_2gev")
# # Slide number 8
# control_plot(list_2gev, "cal_hit_energy", "(200, 0, 100)", "", "cal_energy_2gev")

# # Slide number 9
# control_plot(list_1gev, "tr1_hit_energy", "(200, 0, 5)", "", "tr1_energy_1gev")
# control_plot(list_1gev, "tr2_hit_energy", "(200, 0, 5)", "", "tr2_energy_1gev")
# control_plot(list_1gev, "tr1_hit_y-cal_cluster_y[0]", "(200, -60, 45)", "", "tr1_y_dist_1gev")
# control_plot(list_1gev, "tr2_hit_y-cal_cluster_y[0]", "(200, -60, 45)", "", "tr2_y_dist_1gev")
# # Slide number 10
# control_plot(list_1gev, "cal_hit_energy", "(200, 0, 100)", "", "cal_energy_1gev")


# Slide number 11, 12
# avr_tr_energy(1)
# avr_tr_energy(2)

check_alignment(t_data_1gev)
check_alignment(t_data_2gev)
check_alignment(t_data_3gev)
check_alignment(t_data_4gev)
check_alignment(t_data_5gev)
