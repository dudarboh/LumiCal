from ROOT import TFile, gROOT, TGraphErrors, TH1F
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def residuals_111_before_alignment(input_file, output_file):
    '''
    Mean of residuals before misalignment:
    Tr1: -0.1897
    Tr2: 0.9405
    Cal: -0.7501
    '''
    tree = input_file.data
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


def residuals_111_2points(input_file, output_file):
    tree = input_file.data
    tr1_z = 0.
    tr2_z = 5 * 4.5
    cal_z = 25 * 4.5
    x = np.array([tr1_z, cal_z])
    h_tr2_res = TH1F('h_tr2_res111_2points', 'Tracker2 residuals', 100, -10, 10)

    for event in tree:
        if (event.tr1_n_clusters == 1 and event.tr2_n_clusters == 1 and event.cal_n_clusters == 1
           and 153.1 < event.cal_cluster_rho[0] < 172.3):
            y = np.array([event.tr1_cluster_rho[0], event.cal_cluster_rho[0]])
            x_err = np.array([4.5 / 2, 13.5])
            y_err = np.array([0.9, 0.9])
            track = TGraphErrors(2, x, y, x_err, y_err)
            track.Fit('pol1', "Q")
            fit_func = track.GetFunction('pol1')
            tr2_residual = event.tr2_cluster_rho[0] - fit_func.Eval(tr2_z)
            h_tr2_res.Fill(tr2_residual)

    output_file.cd()
    h_tr2_res.Write()


def residuals_111_cluster(input_file, output_file):
    tree = input_file.data

    h_tr1_res = TH1F('h_tr1_res111_cluster', 'Tracker1 residuals cluster', 200, -20, 20)
    h_tr2_res = TH1F('h_tr2_res111_cluster', 'Tracker2 residuals cluster', 200, -20, 20)

    for event in tree:
        if (event.tr1_n_clusters == 1 and event.tr2_n_clusters == 1 and event.cal_n_clusters == 1
           and 153.1 < event.cal_cluster_rho[0] < 172.3):

            tr1_residual = event.tr1_cluster_rho[0] - event.cal_cluster_rho[0]
            tr2_residual = event.tr2_cluster_rho[0] - event.cal_cluster_rho[0]
            h_tr1_res.Fill(tr1_residual)
            h_tr2_res.Fill(tr2_residual)

    output_file.cd()
    h_tr1_res.Write()
    h_tr2_res.Write()


def plot_backscattered_tracks(input_file):
    tree = input_file.data
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    idx = 0
    for event in tree:
        if (event.tr1_n_clusters == 2 and event.tr2_n_clusters == 2 and event.cal_n_clusters == 1
           and 153.1 < event.cal_cluster_rho[0] < 172.3):
            x = np.linspace(0, 135, 31)

            y0 = event.tr1_cluster_x[1]
            ky = (event.tr2_cluster_x[1] - y0) / 22.5  # mm
            y = y0 + x * ky

            z0 = event.tr1_cluster_y[1]
            kz = (event.tr2_cluster_y[1] - z0) / 22.5  # mm
            z = z0 + x * kz

            ax.plot(x, y, z)
            idx += 1
            if idx == 30:
                break
    ax.set_xlim(0, 135)
    ax.set_ylim(-51, 51)
    ax.set_zlim(153, 173)
    ax.set_xlabel('z (layer), mm')
    ax.set_ylabel('x (sector), mm')
    ax.set_zlabel('y (pad), mm')
    plt.show()


def tr2_efficiency(input_file):
    tree = input_file.data
    eff = 0
    eff_tot = 0
    for event in tree:
        if (event.tr1_n_clusters == 1 and event.cal_n_clusters == 1
           and 153.1 < event.cal_cluster_rho[0] < 172.3
           and -2 < event.tr1_cluster_rho[0] - event.cal_cluster_rho[0] < 2):
            eff_tot += 1
            if event.tr2_n_clusters >= 1 and -2 < event.tr2_cluster_rho[0] - event.cal_cluster_rho[0] < 2:
                eff += 1
            elif event.tr2_n_clusters >= 2 and -2 < event.tr2_cluster_rho[1] - event.cal_cluster_rho[0] < 2:
                print('WATAFAAACK??')

    print('Efficiency of Tracker2:', eff / eff_tot)


def tr1_efficiency(input_file):
    tree = input_file.data
    eff = 0
    eff_tot = 0
    for event in tree:
        if (event.tr2_n_clusters == 1 and event.cal_n_clusters == 1
           and 153.1 < event.cal_cluster_rho[0] < 172.3
           and -2 < event.tr2_cluster_rho[0] - event.cal_cluster_rho[0] < 2):
            eff_tot += 1

            if event.tr1_n_clusters >= 1 and -2 < event.tr1_cluster_rho[0] - event.cal_cluster_rho[0] < 2:
                eff += 1
            elif event.tr1_n_clusters >= 2 and -2 < event.tr1_cluster_rho[1] - event.cal_cluster_rho[0] < 2:
                print('WATAFAAACK??')

    print('Efficiency of Tracker1:', eff / eff_tot)


def merge_pos_change(input_merged_file, input_not_merged_file, output_file):
    tree_merged = input_merged_file.data
    tree_not_merged = input_not_merged_file.data

    h_pos_shift = TH1F('h_pos_shift', 'Merge shift position', 100, -10, 10)

    merged_pos = np.array([])
    not_merged_pos = np.array([])

    print('I am here 1')
    for event1 in tree_merged:
        if event1.cal_n_clusters >= 1:
            merged_pos = np.append(merged_pos, event1.cal_cluster_y[0])

    print('I am here 2')
    for event2 in tree_not_merged:
        if event2.cal_n_clusters >= 1:
            not_merged_pos = np.append(not_merged_pos, event2.cal_cluster_y[0])

    print('I am here 3')
    for idx, pos in enumerate(merged_pos):
        h_pos_shift.Fill(not_merged_pos[idx] - pos)

    output_file.cd()
    h_pos_shift.Write()


def main():

    file_merged = TFile.Open('./extracted_data_merged.root')
    file_not_merged = TFile.Open('./extracted_data_not_merged.root')
    tree = file_merged.data

    output_file = TFile('RENAME.root', 'recreate')

    # plot_backscattered_tracks(file_merged)

    residuals_111_before_alignment(file_merged, output_file)
    # residuals_111_cluster(file_merged, output_file)
    # residuals_111_2points(file_merged, output_file)

    # tr1_efficiency(file_merged)
    # tr2_efficiency(file_merged)

    # merge_pos_change(file_merged, file_not_merged, output_file)

    tree.Draw('cal_cluster_energy[0]>>h_cal_energy', 'cal_n_clusters == 1 && 153.1 < cal_cluster_rho[0] && cal_cluster_rho[0] < 172.3')
    h_cal_energy = gROOT.FindObject('h_cal_energy')
    h_cal_energy.Write()

    tree.Draw('cal_cluster_rho[0]>>h_cal_rho', 'cal_n_clusters == 1 && 153.1 < cal_cluster_rho[0] && cal_cluster_rho[0] < 172.3')
    h_cal_rho = gROOT.FindObject('h_cal_rho')
    h_cal_rho.Write()

    input('Yaay I am finished :3')


gROOT.SetBatch(1)

main()

'''
TO DO LIST:
1)Change alignment!
2) Plot residuals to the 3 point fit
3) Calculate efficiency of Tr2 as distance to
3.1) pol0 fit line through 2 points
3.2) distance to the cluster in calorimeter
4) The same way around the eff for Tr1.
5) Scatter plot of fitted trackes to calorimeter
'''
