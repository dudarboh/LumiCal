from ROOT import TFile, gROOT, TGraphErrors, TH1F, TGraph, TCanvas, TPad, TF1, gStyle, TF1, nullptr, TH2F, TPaletteAxis, gPad, gStyle
import numpy as np

gROOT.SetBatch(0)
gROOT.SetStyle('ATLAS')

# file_data = TFile.Open("../extracted/extracted_data_5gev.root", 'read')
file_data = TFile.Open("../extracted/extracted_data_5gev_merged.root", 'read')
tree_data = file_data.data

output_file = TFile("output.root", "RECREATE")
output_file.cd()


def plot():
    h1 = TH1F('h1', 'Clst1', 8, 0, 8)
    h2 = TH1F('h2', 'Clst2', 8, 0, 8)
    # for event in tree_data:
    tree_data.Draw("tr1_n_hits>>h1", "", "histo")
    h1.SetLineColor(1)

    tree_data.Draw("tr2_n_hits>>h2", "", "histosame")
    h2.SetLineColor(2)

    input("wait")


def identification():
    x = np.arange(0, 10, 0.1)
    y1 = []
    y2 = []

    tree_data.Draw(">>n_events1", "cal_n_clusters>0")
    n_events1 = gROOT.FindObject("n_events1")
    n_ev1 = n_events1.GetN()

    tree_data.Draw(">>n_events2", "cal_n_clusters>1")
    n_events2 = gROOT.FindObject("n_events2")
    n_ev2 = n_events2.GetN()


    for acc in x:
        print("Obtaining {} point".format(acc))
        tree_data.Draw(">>list1", "(Sum$(abs(cal_cluster_y[0]-tr1_hit_y)<{0}) > 0) && (Sum$(abs(cal_cluster_y[0]-tr2_hit_y)<{0}) > 0)".format(acc))
        # tree_data.Draw(">>list2", "Sum$(abs(cal_cluster_y[0]-tr2_hit_y)<{}) > 0".format(acc))

        l1 = gROOT.FindObject("list1")
        # l2 = gROOT.FindObject("list2")

        y1.append(l1.GetN() / n_ev1 * 100.)
        # y2.append(l2.GetN() / n_ev1 * 100.)

    gr1 = TGraphErrors(len(x), x, np.array(y1), nullptr, nullptr)
    gr1.Draw("APL")
    gr1.SetTitle("Hits in both trackers")
    gr1.GetXaxis().SetTitle("Accepted track distance, mm")
    gr1.GetYaxis().SetTitle("Accepted events, %")

    # gr2 = TGraphErrors(len(x), x, np.array(y2), nullptr, nullptr)
    # gr2.Draw("PLsame")
    # gr2.SetTitle("Tracker 2 hits")
    # gr2.SetMarkerColor(2)
    # gr2.SetLineColor(2)

    input("wait")


def method1(fit_func="pol0"):
    """Fit the horizontal line through reconstructed shower and hit in tracker2.
    Calculate expected position of hit in tracker1 from the fit of those 2 points.
    Plot histo of residuals of expected and measured positions.
    Plot efficiency and purity points for various cut-off values on the residual."""
    h_reco_gen = TH1F("h_reco_gen", "title", 2000, -10., 10.)
    h_reco = TH1F("h_reco", "title", 2000, -10., 10.)
    n_gen = 0

    n_events = tree_data.GetEntries()
    x = np.array([6., 24.])
    y_err = np.array([0.52, 0.44])

    for i, event in enumerate(tree_data):
        if i % 1000 == 0:
            if i == 100000:
                break
            print(i, "event out of", n_events)

        # Good event with clear electron generated and detected by calorimeter and trackers
        if not (event.cal_n_clusters > 0
                and event.cal_cluster_energy[0] > 150
                and 130 < event.cal_cluster_y[0] < 150
                and event.tr1_n_hits > 0
                and event.tr2_n_hits > 0):
            continue

        n_gen += 1

        # Do the algorithm for generated track
        # This is algorithm
        y = np.array([event.tr2_hit_y[0], event.cal_cluster_y[0]])
        gr_track = TGraphErrors(2, x, y, nullptr, y_err)
        gr_track.Fit(fit_func, "Q")
        residual = event.tr1_hit_y[0] - gr_track.GetFunction(fit_func).Eval(1.)
        h_reco_gen.Fill(residual)

        for y1 in event.tr1_hit_y:
            for y2 in event.tr2_hit_y:
                for y3 in event.cal_cluster_y:
                    # This is algorithm
                    y = np.array([y2, y3])
                    gr_track = TGraphErrors(2, x, y, nullptr, y_err)
                    gr_track.Fit(fit_func, "Q")
                    residual = y1 - gr_track.GetFunction(fit_func).Eval(1.)

                    h_reco.Fill(residual)

    h_reco_gen.Write("method1_reco_gen_" + fit_func)
    h_reco.Write("method1_reco_" + fit_func)
    # Get the graph
    cut_off = np.arange(0.01, 10.01, 0.01)
    eff = []
    purity = []

    for cut in cut_off:
        bmin = h_reco_gen.GetXaxis().FindBin(0.548 - cut)
        bmax = h_reco_gen.GetXaxis().FindBin(0.548 + cut)
        n_reco_gen = h_reco_gen.Integral(bmin, bmax)
        n_reco = h_reco.Integral(bmin, bmax)

        if n_gen != 0 and n_reco != 0:
            eff.append(n_reco_gen / n_gen * 100)
            purity.append(n_reco_gen / n_reco * 100)

    gr = TGraph(len(cut_off), np.array(eff), np.array(purity))
    gr.GetXaxis().SetTitle("Efficiency, %")
    gr.GetYaxis().SetTitle("Purity, %")
    gr.SetTitle("Tr2 Cal fit " + fit_func)

    gr.Write("method1_plot_" + fit_func)

    return gr


def method2(fit_func="pol0"):
    """Fit the horizontal line through reconstructed shower and hit in tracker2.
    Calculate expected position of hit in tracker1 from the fit of those 2 points.
    Plot histo of residuals of expected and measured positions.
    Plot efficiency and purity points for various cut-off values on the residual."""
    h_reco_gen = TH1F("h_reco_gen", "title", 2000, -10., 10.)
    h_reco = TH1F("h_reco", "title", 2000, -10., 10.)
    n_gen = 0

    n_events = tree_data.GetEntries()
    x = np.array([1., 24.])
    y_err = np.array([0.52, 0.44])

    for i, event in enumerate(tree_data):
        if i % 1000 == 0:
            if i == 100000:
                break
            print(i, "event out of", n_events)

        # Good event with clear electron generated and detected by calorimeter and trackers
        if not (event.cal_n_clusters > 0
                and event.cal_cluster_energy[0] > 150
                and 130 < event.cal_cluster_y[0] < 150
                and event.tr1_n_hits > 0
                and event.tr2_n_hits > 0):
            continue

        n_gen += 1

        # Do the algorithm for generated track
        # This is algorithm
        y = np.array([event.tr1_hit_y[0], event.cal_cluster_y[0]])
        gr_track = TGraphErrors(2, x, y, nullptr, y_err)
        gr_track.Fit(fit_func, "Q")
        residual = event.tr2_hit_y[0] - gr_track.GetFunction(fit_func).Eval(6.)
        h_reco_gen.Fill(residual)

        for y1 in event.tr1_hit_y:
            for y2 in event.tr2_hit_y:
                for y3 in event.cal_cluster_y:
                    # This is algorithm
                    y = np.array([y1, y3])
                    gr_track = TGraphErrors(2, x, y, nullptr, y_err)
                    gr_track.Fit(fit_func, "Q")
                    residual = y2 - gr_track.GetFunction(fit_func).Eval(6.)

                    h_reco.Fill(residual)

    h_reco_gen.Write("method2_reco_gen_" + fit_func)
    h_reco.Write("method2_reco_" + fit_func)
    # Get the graph
    cut_off = np.arange(0.01, 10.01, 0.01)
    eff = []
    purity = []

    for cut in cut_off:
        bmin = h_reco_gen.GetXaxis().FindBin(0.254 - cut)
        bmax = h_reco_gen.GetXaxis().FindBin(0.254 + cut)
        n_reco_gen = h_reco_gen.Integral(bmin, bmax)
        n_reco = h_reco.Integral(bmin, bmax)

        if n_gen != 0 and n_reco != 0:
            eff.append(n_reco_gen / n_gen * 100)
            purity.append(n_reco_gen / n_reco * 100)

    gr = TGraph(len(cut_off), np.array(eff), np.array(purity))
    gr.GetXaxis().SetTitle("Efficiency, %")
    gr.GetYaxis().SetTitle("Purity, %")
    gr.SetTitle("Tr1 Cal fit " + fit_func)

    gr.Write("method2_plot_" + fit_func)

    return gr


def method3(fit_func="pol0"):
    """Fit the horizontal line through hits in tracker1 and tracker2.
    Calculate expected position of hit in the calorimeter from the fit of those 2 points.
    Plot histo of residuals of expected and measured positions of those hits.
    Plot efficiency and purity points for various cut-off values on the residual."""
    h_reco_gen = TH1F("h_reco_gen", "title", 2000, -10., 10.)
    h_reco = TH1F("h_reco", "title", 2000, -10., 10.)
    n_gen = 0

    n_events = tree_data.GetEntries()
    x = np.array([1., 6.])
    y_err = np.array([0.52, 0.52])

    for i, event in enumerate(tree_data):
        if i % 1000 == 0:
            if i == 100000:
                break
            print(i, "event out of", n_events)

        # Good event with clear electron generated and detected by calorimeter and trackers
        if not (event.cal_n_clusters > 0
                and event.cal_cluster_energy[0] > 150
                and 130 < event.cal_cluster_y[0] < 150
                and event.tr1_n_hits > 0
                and event.tr2_n_hits > 0):
            continue

        n_gen += 1

        # Do the algorithm for generated track
        # This is algorithm
        y = np.array([event.tr1_hit_y[0], event.tr2_hit_y[0]])
        gr_track = TGraphErrors(2, x, y, nullptr, y_err)
        gr_track.Fit(fit_func, "Q")
        residual = event.cal_cluster_y[0] - gr_track.GetFunction(fit_func).Eval(24.)
        h_reco_gen.Fill(residual)

        for y1 in event.tr1_hit_y:
            for y2 in event.tr2_hit_y:
                for y3 in event.cal_cluster_y:
                    # This is algorithm
                    y = np.array([y1, y2])
                    gr_track = TGraphErrors(2, x, y, nullptr, y_err)
                    gr_track.Fit(fit_func, "Q")
                    residual = y3 - gr_track.GetFunction(fit_func).Eval(24.)

                    h_reco.Fill(residual)

    h_reco_gen.Write("method3_reco_gen_" + fit_func)
    h_reco.Write("method3_reco_" + fit_func)
    # Get the graph
    cut_off = np.arange(0.01, 10.01, 0.01)
    eff = []
    purity = []

    mean_reco_gen = h_reco_gen.GetMean()
    mean_reco = h_reco.GetMean()
    for cut in cut_off:
        bmin_reco_gen = h_reco_gen.GetXaxis().FindBin(mean_reco_gen - cut)
        bmax_reco_gen = h_reco_gen.GetXaxis().FindBin(mean_reco_gen + cut)
        n_reco_gen = h_reco_gen.Integral(bmin_reco_gen, bmax_reco_gen)

        bmin_reco = h_reco.GetXaxis().FindBin(mean_reco - cut)
        bmax_reco = h_reco.GetXaxis().FindBin(mean_reco + cut)
        n_reco = h_reco.Integral(bmin_reco, bmax_reco)

        if n_gen != 0 and n_reco != 0:
            eff.append(n_reco_gen / n_gen * 100)
            purity.append(n_reco_gen / n_reco * 100)

    gr = TGraph(len(cut_off), np.array(eff), np.array(purity))
    gr.GetXaxis().SetTitle("Efficiency, %")
    gr.GetYaxis().SetTitle("Purity, %")
    gr.SetTitle("Tr1 Tr2 fit " + fit_func)

    gr.Write("method3_plot_" + fit_func)

    return gr


def method4(fit_func="pol0"):
    """Fit the horizontal line through hits in tracker1, tracker2 and calorimeter.
    Calculate expected position of hit in the calorimeter from the fit of those 3 points.
    Plot histo of residuals of expected and measured positions of the shower.
    Plot efficiency and purity points for various cut-off values on the residual."""
    h_reco_gen = TH1F("h_reco_gen", "title", 2000, -10., 10.)
    h_reco = TH1F("h_reco", "title", 2000, -10., 10.)
    n_gen = 0

    n_events = tree_data.GetEntries()
    x = np.array([1., 6., 24.])
    y_err = np.array([0.52, 0.52, 0.44])

    for i, event in enumerate(tree_data):
        if i % 1000 == 0:
            if i == 100000:
                break
            print(i, "event out of", n_events)

        # Good event with clear electron generated and detected by calorimeter and trackers
        if not (event.cal_n_clusters > 0
                and event.cal_cluster_energy[0] > 150
                and 130 < event.cal_cluster_y[0] < 150
                and event.tr1_n_hits > 0
                and event.tr2_n_hits > 0):
            continue

        n_gen += 1

        # Do the algorithm for generated track
        # This is algorithm
        y = np.array([event.tr1_hit_y[0], event.tr2_hit_y[0], event.cal_cluster_y[0]])
        gr_track = TGraphErrors(3, x, y, nullptr, y_err)
        gr_track.Fit(fit_func, "Q")
        residual = event.cal_cluster_y[0] - gr_track.GetFunction(fit_func).Eval(24.)
        h_reco_gen.Fill(residual)

        for y1 in event.tr1_hit_y:
            for y2 in event.tr2_hit_y:
                for y3 in event.cal_cluster_y:
                    # This is algorithm
                    y = np.array([y1, y2, y3])
                    gr_track = TGraphErrors(3, x, y, nullptr, y_err)
                    gr_track.Fit(fit_func, "Q")
                    residual = y3 - gr_track.GetFunction(fit_func).Eval(24.)

                    h_reco.Fill(residual)

    h_reco_gen.Write("method4_reco_gen_" + fit_func)
    h_reco.Write("method4_reco_" + fit_func)
    # Get the graph
    cut_off = np.arange(0.01, 10.01, 0.01)
    eff = []
    purity = []

    mean_reco_gen = h_reco_gen.GetMean()
    mean_reco = h_reco.GetMean()
    for cut in cut_off:
        bmin_reco_gen = h_reco_gen.GetXaxis().FindBin(mean_reco_gen - cut)
        bmax_reco_gen = h_reco_gen.GetXaxis().FindBin(mean_reco_gen + cut)
        n_reco_gen = h_reco_gen.Integral(bmin_reco_gen, bmax_reco_gen)

        bmin_reco = h_reco.GetXaxis().FindBin(mean_reco - cut)
        bmax_reco = h_reco.GetXaxis().FindBin(mean_reco + cut)
        n_reco = h_reco.Integral(bmin_reco, bmax_reco)

        if n_gen != 0 and n_reco != 0:
            eff.append(n_reco_gen / n_gen * 100)
            purity.append(n_reco_gen / n_reco * 100)

    gr = TGraph(len(cut_off), np.array(eff), np.array(purity))
    gr.GetXaxis().SetTitle("Efficiency, %")
    gr.GetYaxis().SetTitle("Purity, %")
    gr.SetTitle("Tr1 Tr2 Cal fit " + fit_func)

    gr.Write("method4_plot_" + fit_func)

    return gr


gr1 = method1()
gr1.Draw("AP")

gr2 = method2()
gr2.SetMarkerColor(2)
gr2.SetMarkerStyle(21)
gr2.Draw("Psame")

gr3 = method3()
gr3.SetMarkerColor(3)
gr3.SetMarkerStyle(22)
gr3.Draw("Psame")

gr4 = method4()
gr4.SetMarkerColor(4)
gr4.SetMarkerStyle(23)
gr4.Draw("Psame")
###

# gr5 = method1("pol1")
# gr5.SetMarkerColor(5)
# gr5.SetMarkerStyle(24)
# gr5.Draw("Psame")

# gr6 = method2("pol1")
# gr6.SetMarkerColor(6)
# gr6.SetMarkerStyle(25)
# gr6.Draw("Psame")

# gr7 = method3("pol1")
# gr7.SetMarkerColor(7)
# gr7.SetMarkerStyle(26)
# gr7.Draw("Psame")

# gr8 = method4("pol1")
# gr8.SetMarkerColor(8)
# gr8.SetMarkerStyle(27)
# gr8.Draw("Psame")


input("wait")
