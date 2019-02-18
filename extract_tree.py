from ROOT import TFile, TTree
import array
import time

from signal_readout import Signal, position, calib_energy
from clustering import Cluster, set_neighbors, set_clusters, merge_clusters


def extract_signal(event):
    id_arr = event.apv_id
    channel_arr = event.apv_ch
    signal_arr = event.apv_signal_maxfit
    apv_nn_output = event.apv_nn_output
    apv_fit_tau = event.apv_fit_tau
    apv_fit_t0 = event.apv_fit_t0
    apv_bint1 = event.apv_bint1

    signals_calorimeter = []
    signals_tracker1 = []
    signals_tracker2 = []

    for hit in range(len(id_arr)):
        if (apv_fit_tau[hit] < 1 or apv_fit_tau[hit] > 3
           or signal_arr[hit] > 2000.
           or apv_fit_t0[hit] < (apv_bint1[hit] - 2.7)
           or apv_fit_t0[hit] > (apv_bint1[hit] - 0.5)):
            continue

        sector, pad, layer = position(id_arr[hit], channel_arr[hit])
        energy = calib_energy(id_arr[hit], signal_arr[hit])

        if (sector == 0 or sector == 3 or layer == 7
           or (sector == 1 and pad < 20)
           or (sector == 2 and pad < 20)
           or sector < 0  # This one is changed due to python C++ difference in %.
           or (layer < 2 and (signal_arr[hit] < 0. or apv_nn_output[hit] < 0.5))
           or (layer >= 2 and (energy < 1.4 or apv_nn_output[hit] < 0.5))):
            continue

        if layer == 0:
            data_list = signals_tracker1
        elif layer == 1:
            data_list = signals_tracker2
        else:
            data_list = signals_calorimeter

        for existed in data_list:
            if (existed.sector, existed.pad) == (sector, pad):
                existed.energy += energy
                break
        else:
            data_list.append(Signal(sector, pad, layer, energy))

    signals_calorimeter = sorted(signals_calorimeter, key=lambda x: x.energy, reverse=True)
    signals_tracker1 = sorted(signals_tracker1, key=lambda x: x.energy, reverse=True)
    signals_tracker2 = sorted(signals_tracker2, key=lambda x: x.energy, reverse=True)

    return signals_tracker1, signals_tracker2, signals_calorimeter


def clustering_in_towers(signals_list, merge, det):
    clusters_list = []

    if len(signals_list) == 0:
        return clusters_list

    set_neighbors(signals_list)
    set_clusters(signals_list)

    n_clusters = max([signal.cluster for signal in signals_list]) + 1
    for cluster in range(n_clusters):
        if det == 'Cal':
            clusters_list.append(Cluster(signals_list, cluster, 'logW'))
        else:
            clusters_list.append(Cluster(signals_list, cluster, 'Energy'))

    if merge == 'on':
        merge_clusters(signals_list, clusters_list)

    clusters_list = sorted(clusters_list, key=lambda x: x.energy, reverse=True)

    return clusters_list


def main():

    start_time = time.time()

    file = TFile.Open('./trees/run741_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root')
    tree = file.apv_reco
    n_events = tree.GetEntries()

    output_file = TFile('extracted_data_RENAME.root', 'recreate')
    output_tree = TTree('data', 'Extracted data')

    tr1_n_signals = array.array('i', [0])
    tr1_signal_pad = array.array('i', [0] * 256)
    tr1_signal_sector = array.array('i', [0] * 256)
    tr1_signal_rho = array.array('f', [0.0] * 256)
    tr1_signal_x = array.array('f', [0.0] * 256)
    tr1_signal_y = array.array('f', [0.0] * 256)
    tr1_signal_energy = array.array('f', [0.0] * 256)
    tr1_n_clusters = array.array('i', [0])
    tr1_cluster_pad = array.array('f', [0.0] * 256)
    tr1_cluster_sector = array.array('f', [0.0] * 256)
    tr1_cluster_rho = array.array('f', [0.0] * 256)
    tr1_cluster_x = array.array('f', [0.0] * 256)
    tr1_cluster_y = array.array('f', [0.0] * 256)
    tr1_cluster_energy = array.array('f', [0.0] * 256)
    tr1_cluster_n_pads = array.array('f', [0.0] * 256)

    tr2_n_signals = array.array('i', [0])
    tr2_signal_pad = array.array('i', [0] * 256)
    tr2_signal_sector = array.array('i', [0] * 256)
    tr2_signal_rho = array.array('f', [0.0] * 256)
    tr2_signal_x = array.array('f', [0.0] * 256)
    tr2_signal_y = array.array('f', [0.0] * 256)
    tr2_signal_energy = array.array('f', [0.0] * 256)
    tr2_n_clusters = array.array('i', [0])
    tr2_cluster_pad = array.array('f', [0.0] * 256)
    tr2_cluster_sector = array.array('f', [0.0] * 256)
    tr2_cluster_rho = array.array('f', [0.0] * 256)
    tr2_cluster_x = array.array('f', [0.0] * 256)
    tr2_cluster_y = array.array('f', [0.0] * 256)
    tr2_cluster_energy = array.array('f', [0.0] * 256)
    tr2_cluster_n_pads = array.array('f', [0.0] * 256)

    cal_n_signals = array.array('i', [0])
    cal_signal_pad = array.array('i', [0] * 256)
    cal_signal_sector = array.array('i', [0] * 256)
    cal_signal_rho = array.array('f', [0.0] * 256)
    cal_signal_x = array.array('f', [0.0] * 256)
    cal_signal_y = array.array('f', [0.0] * 256)
    cal_signal_energy = array.array('f', [0.0] * 256)
    cal_n_clusters = array.array('i', [0])
    cal_cluster_pad = array.array('f', [0.0] * 256)
    cal_cluster_sector = array.array('f', [0.0] * 256)
    cal_cluster_rho = array.array('f', [0.0] * 256)
    cal_cluster_x = array.array('f', [0.0] * 256)
    cal_cluster_y = array.array('f', [0.0] * 256)
    cal_cluster_energy = array.array('f', [0.0] * 256)
    cal_cluster_n_pads = array.array('f', [0.0] * 256)

    output_tree.Branch('tr1_n_signals', tr1_n_signals, 'tr1_n_signals/I')
    output_tree.Branch('tr1_signal_pad', tr1_signal_pad, 'tr1_signal_pad[tr1_n_signals]/I')
    output_tree.Branch('tr1_signal_sector', tr1_signal_sector, 'tr1_signal_sector[tr1_n_signals]/I')
    output_tree.Branch('tr1_signal_rho', tr1_signal_rho, 'tr1_signal_rho[tr1_n_signals]/F')
    output_tree.Branch('tr1_signal_x', tr1_signal_x, 'tr1_signal_x[tr1_n_signals]/F')
    output_tree.Branch('tr1_signal_y', tr1_signal_y, 'tr1_signal_y[tr1_n_signals]/F')
    output_tree.Branch('tr1_signal_energy', tr1_signal_energy, 'tr1_signal_energy[tr1_n_signals]/F')
    output_tree.Branch('tr1_n_clusters', tr1_n_clusters, 'tr1_n_clusters/I')
    output_tree.Branch('tr1_cluster_pad', tr1_cluster_pad, 'tr1_cluster_pad[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_sector', tr1_cluster_sector, 'tr1_cluster_sector[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_rho', tr1_cluster_rho, 'tr1_cluster_rho[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_x', tr1_cluster_x, 'tr1_cluster_x[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_y', tr1_cluster_y, 'tr1_cluster_y[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_energy', tr1_cluster_energy, 'tr1_cluster_energy[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_n_pads', tr1_cluster_n_pads, 'tr1_cluster_n_pads[tr1_n_clusters]/F')

    output_tree.Branch('tr2_n_signals', tr2_n_signals, 'tr2_n_signals/I')
    output_tree.Branch('tr2_signal_pad', tr2_signal_pad, 'tr2_signal_pad[tr2_n_signals]/I')
    output_tree.Branch('tr2_signal_sector', tr2_signal_sector, 'tr2_signal_sector[tr2_n_signals]/I')
    output_tree.Branch('tr2_signal_rho', tr2_signal_rho, 'tr2_signal_rho[tr2_n_signals]/F')
    output_tree.Branch('tr2_signal_x', tr2_signal_x, 'tr2_signal_x[tr2_n_signals]/F')
    output_tree.Branch('tr2_signal_y', tr2_signal_y, 'tr2_signal_y[tr2_n_signals]/F')
    output_tree.Branch('tr2_signal_energy', tr2_signal_energy, 'tr2_signal_energy[tr2_n_signals]/F')
    output_tree.Branch('tr2_n_clusters', tr2_n_clusters, 'tr2_n_clusters/I')
    output_tree.Branch('tr2_cluster_pad', tr2_cluster_pad, 'tr2_cluster_pad[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_sector', tr2_cluster_sector, 'tr2_cluster_sector[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_rho', tr2_cluster_rho, 'tr2_cluster_rho[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_x', tr2_cluster_x, 'tr2_cluster_x[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_y', tr2_cluster_y, 'tr2_cluster_y[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_energy', tr2_cluster_energy, 'tr2_cluster_energy[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_n_pads', tr2_cluster_n_pads, 'tr2_cluster_n_pads[tr2_n_clusters]/F')

    output_tree.Branch('cal_n_signals', cal_n_signals, 'cal_n_signals/I')
    output_tree.Branch('cal_signal_pad', cal_signal_pad, 'cal_signal_pad[cal_n_signals]/I')
    output_tree.Branch('cal_signal_sector', cal_signal_sector, 'cal_signal_sector[cal_n_signals]/I')
    output_tree.Branch('cal_signal_rho', cal_signal_rho, 'cal_signal_rho[cal_n_signals]/F')
    output_tree.Branch('cal_signal_x', cal_signal_x, 'cal_signal_x[cal_n_signals]/F')
    output_tree.Branch('cal_signal_y', cal_signal_y, 'cal_signal_y[cal_n_signals]/F')
    output_tree.Branch('cal_signal_energy', cal_signal_energy, 'cal_signal_energy[cal_n_signals]/F')
    output_tree.Branch('cal_n_clusters', cal_n_clusters, 'cal_n_clusters/I')
    output_tree.Branch('cal_cluster_pad', cal_cluster_pad, 'cal_cluster_pad[cal_n_clusters]/F')
    output_tree.Branch('cal_cluster_sector', cal_cluster_sector, 'cal_cluster_sector[cal_n_clusters]/F')
    output_tree.Branch('cal_cluster_rho', cal_cluster_rho, 'cal_cluster_rho[cal_n_clusters]/F')
    output_tree.Branch('cal_cluster_x', cal_cluster_x, 'cal_cluster_x[cal_n_clusters]/F')
    output_tree.Branch('cal_cluster_y', cal_cluster_y, 'cal_cluster_y[cal_n_clusters]/F')
    output_tree.Branch('cal_cluster_energy', cal_cluster_energy, 'cal_cluster_energy[cal_n_clusters]/F')
    output_tree.Branch('cal_cluster_n_pads', cal_cluster_n_pads, 'cal_cluster_n_pads[cal_n_clusters]/F')

    for idx, event in enumerate(tree):
        # if idx == 3000:
        #    break

        if idx % (1000) == 0:
            time_min = (time.time() - start_time) // 60
            time_sec = (time.time() - start_time) % 60
            print('%d/%d events' % (idx, n_events), end=' ')
            print('%d min' % time_min, end=' ')
            print('%d sec' % time_sec)

        # Extract data.
        signals_tr1, signals_tr2, signals_cal = extract_signal(event)

        tr1_n_signals[0] = len(signals_tr1)
        for i, signal in enumerate(signals_tr1):
            tr1_signal_pad[i] = signal.pad
            tr1_signal_sector[i] = signal.sector
            tr1_signal_rho[i] = signal.rho
            tr1_signal_x[i] = signal.x
            tr1_signal_y[i] = signal.y
            tr1_signal_energy[i] = signal.energy

        tr2_n_signals[0] = len(signals_tr2)
        for i, signal in enumerate(signals_tr2):
            tr2_signal_pad[i] = signal.pad
            tr2_signal_sector[i] = signal.sector
            tr2_signal_rho[i] = signal.rho
            tr2_signal_x[i] = signal.x
            tr2_signal_y[i] = signal.y
            tr2_signal_energy[i] = signal.energy

        cal_n_signals[0] = len(signals_cal)
        for i, signal in enumerate(signals_cal):
            cal_signal_pad[i] = signal.pad
            cal_signal_sector[i] = signal.sector
            cal_signal_rho[i] = signal.rho
            cal_signal_x[i] = signal.x
            cal_signal_y[i] = signal.y
            cal_signal_energy[i] = signal.energy

        clusters_cal = clustering_in_towers(signals_cal, merge='on', det='Cal')
        clusters_tr1 = clustering_in_towers(signals_tr1, merge='off', det='Tr1')
        clusters_tr2 = clustering_in_towers(signals_tr2, merge='off', det='Tr2')
        if len(clusters_cal) != 0:
            clusters_tr1 = sorted(clusters_tr1, key=lambda x: abs(x.rho - clusters_cal[0].rho))
            clusters_tr2 = sorted(clusters_tr2, key=lambda x: abs(x.rho - clusters_cal[0].rho))

        tr1_n_clusters[0] = len(clusters_tr1)
        for i, cluster in enumerate(clusters_tr1):
            tr1_cluster_pad[i] = cluster.pad
            tr1_cluster_sector[i] = cluster.sector
            tr1_cluster_rho[i] = cluster.rho
            tr1_cluster_x[i] = cluster.x
            tr1_cluster_y[i] = cluster.y
            tr1_cluster_energy[i] = cluster.energy
            tr1_cluster_n_pads[i] = cluster.n_pads

        tr2_n_clusters[0] = len(clusters_tr2)
        for i, cluster in enumerate(clusters_tr2):
            tr2_cluster_pad[i] = cluster.pad
            tr2_cluster_sector[i] = cluster.sector
            tr2_cluster_rho[i] = cluster.rho
            tr2_cluster_x[i] = cluster.x
            tr2_cluster_y[i] = cluster.y
            tr2_cluster_energy[i] = cluster.energy
            tr2_cluster_n_pads[i] = cluster.n_pads

        cal_n_clusters[0] = len(clusters_cal)
        for i, cluster in enumerate(clusters_cal):
            cal_cluster_pad[i] = cluster.pad
            cal_cluster_sector[i] = cluster.sector
            cal_cluster_rho[i] = cluster.rho
            cal_cluster_x[i] = cluster.x
            cal_cluster_y[i] = cluster.y
            cal_cluster_energy[i] = cluster.energy
            cal_cluster_n_pads[i] = cluster.n_pads

        output_tree.Fill()
    output_tree.Write()
    output_file.Close()


main()
