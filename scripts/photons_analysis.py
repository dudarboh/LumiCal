from ROOT import TFile, gROOT, TGraphErrors, TH1F, TGraph, TCanvas, TPad, TF1, gStyle, TF1, nullptr, TH2F, TPaletteAxis, gPad, gStyle
import numpy as np

gROOT.SetBatch(0)
gROOT.SetStyle('ATLAS')

# file_data = TFile.Open("../extracted_root_files/extracted_data_photons_5gev.root", 'read')
file_data = TFile.Open("../extracted_root_files/extracted_data_5gev.root", 'read')
# file_data = TFile.Open("../extracted_root_files/extracted_data_photons_eveto_5gev.root", 'read')
tree_data = file_data.data

output_file = TFile.Open("./results_photons.root", "RECREATE")
output_file.cd()


def plot():
    gStyle.SetPalette(1)
    tree_data.Draw("tr1_cluster_energy>>h(200,0,10)", "", "histo")
    tree_data.Draw("tr1_hit_energy>>h1(200,0,10)", "", "histosame")

    h = gROOT.FindObject("h")
    h.DrawNormalized()

    h1 = gROOT.FindObject("h1")
    h1.DrawNormalized("same")
    input("wait")


def main():
    histo = TH2F('histo', 'histo', 4, 0, 4, 64, 0, 64)
    histo.SetContour(50)

    canvas = TCanvas('name', 'title', 1024, 768)
    canvas.cd()
    counter = 0
    for i, event in enumerate(tree_data):
        if counter == 20:
            break
        energy_in_cal = sum([en for en in event.cal_hit_energy])
        if energy_in_cal < 150:
            continue
        print(energy_in_cal)

        tree_data.Draw("cal_hit_pad:cal_hit_sector>>histo", "cal_hit_energy*(Entry$ == {})".format(i), "colz")
        canvas.Update()
        palette = histo.GetListOfFunctions().FindObject("palette")
        palette.SetX1NDC(0.88)
        palette.SetX2NDC(0.93)
        palette.SetY1NDC(0.2)
        palette.SetY2NDC(0.9)

        histo.SetMinimum(0.)
        histo.SetMaximum(200.)

        canvas.Print("./events_pics/{}".format(i) + '.png')
        counter += 1

    # h1 = TH1F('h1', 'h1', 64, 80., 195.2)
    # h2 = TH1F('h2', 'h2', 64, 80., 195.2)
    # h3 = TH1F('h3', 'h3', 64, 80., 195.2)
    # h4 = TH1F('h4', 'h4', 64, 80., 195.2)
    # # for event in tree_data:
    # tree_data.Draw("cal_cluster_y[0]>>h1", "cal_n_hits>0", "histo")
    # h1.SetLineColor(1)

    # tree_data.Draw("cal_cluster_y[1]>>h2", "cal_n_hits>0", "histosame")
    # h2.SetLineColor(2)

    # tree_data.Draw("cal_cluster_y[2]>>h3", "cal_n_hits>0", "histosame")
    # h3.SetLineColor(3)

    # tree_data.Draw("cal_cluster_y[3]>>h4", "cal_n_hits>0", "histosame")
    # h4.SetLineColor(4)

    # input("wait")


# main()
plot()
