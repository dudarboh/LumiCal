from ROOT import TFile, TTree
import array
import time
import math
import random
from store_file import Hit, HitMC, set_towers, set_clusters, bad_pad


def extract_hits(event):
    id_arr = event.apv_id
    channel_arr = event.apv_ch
    signal_arr = event.apv_signal_maxfit
    apv_nn_output = event.apv_nn_output
    apv_fit_tau = event.apv_fit_tau
    apv_fit_t0 = event.apv_fit_t0
    apv_bint1 = event.apv_bint1

    hits_calorimeter = []
    hits_tracker1 = []
    hits_tracker2 = []

    n_hits = len(id_arr)
    for i in range(n_hits):
        if (apv_fit_tau[i] < 1 or apv_fit_tau[i] > 3
           or signal_arr[i] > 2000.
           or apv_fit_t0[i] < (apv_bint1[i] - 2.7)
           or apv_fit_t0[i] > (apv_bint1[i] - 0.5)):
            continue

        hit = Hit(id_arr[i], channel_arr[i], signal_arr[i])

        if (hit.sector == 0
           or hit.sector == 3
           or hit.layer == 7
           or hit.pad < 0
           or hit.sector < 0  # This one is changed due to python C++ difference in %
           or (hit.sector == 1 and hit.pad < 20)
           or (hit.sector == 2 and hit.pad < 20)
           or (hit.layer < 2 and (signal_arr[i] < 0. or apv_nn_output[i] < 0.5))
           or (hit.layer >= 2 and (hit.energy < 1.4 or apv_nn_output[i] < 0.5))
           or bad_pad(hit.sector, hit.pad, hit.layer)):
            continue

        if hit.layer == 0:
            hits_tracker1.append(hit)
        elif hit.layer == 1:
            hits_tracker2.append(hit)
        else:
            hits_calorimeter.append(hit)

    return hits_tracker1, hits_tracker2, hits_calorimeter


def extract_mc(event):

    vx = event.vX
    vy = event.vY
    px = event.GetLeaf("Tracks.pX").GetValue()
    py = event.GetLeaf("Tracks.pY").GetValue()
    pz = event.GetLeaf("Tracks.pZ").GetValue()
    # Calculated as averaged through all hits at certain zeds
    z_tr1 = 3300.511
    z_tr2 = 3325.513
    z_cal = 3384.033

    # +177.2 is to convert y from montecarto to my coordinate system
    x_tr1 = vx + px * z_tr1 / pz
    y_tr1 = vy + py * z_tr1 / pz + 177.2
    rho_tr1 = (x_tr1**2 + y_tr1**2)**0.5

    x_tr2 = vx + px * z_tr2 / pz
    y_tr2 = vy + py * z_tr2 / pz + 177.2
    rho_tr2 = (x_tr2**2 + y_tr2**2)**0.5

    x_cal = vx + px * z_cal / pz
    y_cal = vy + py * z_cal / pz + 177.2
    rho_cal = (x_cal**2 + y_cal**2)**0.5

    true_hits = [(x_tr1, y_tr1, rho_tr1), (x_tr2, y_tr2, rho_tr2), (x_cal, y_cal, rho_cal)]

    n_hits = event.numHits
    hits_calorimeter = []
    hits_tracker1 = []
    hits_tracker2 = []

    for i in range(n_hits):

        cell_id = event.GetLeaf("Hits.cellID").GetValue(i)
        energy_in_mev = event.GetLeaf("Hits.eHit").GetValue(i)

        hit = HitMC(cell_id, energy_in_mev)

        # Calorimeter efficiency simulation
        S0 = 0.819
        p1 = 2.166
        p0 = 0.999 / 2.
        if random.random() > (1 + math.erf((hit.energy - S0) / p1)) * p0 and hit.layer > 1:
            continue

        if (hit.sector == 0
           or hit.sector == 3
           or hit.layer == 7
           or hit.pad < 0
           or hit.sector < 0
           or (hit.sector == 1 and hit.pad < 20)
           or (hit.sector == 2 and hit.pad < 20)
           or (hit.layer >= 2 and hit.energy < 1.4)
           or bad_pad(hit.sector, hit.pad, hit.layer)):
            continue

        if hit.layer == 0:
            hits_tracker1.append(hit)
        elif hit.layer == 1:
            hits_tracker2.append(hit)
        else:
            hits_calorimeter.append(hit)

    return hits_tracker1, hits_tracker2, hits_calorimeter, true_hits


def align_detector(hits_tr1, hits_tr2, hits_cal):
    tr1_shift = -0.16287581540422025
    tr2_shift = 0.9707477456103106
    cal_shift = -0.807871930206062
    for hit in hits_tr1:
        hit.y -= tr1_shift
    for hit in hits_tr2:
        hit.y -= tr2_shift
    for hit in hits_cal:
        hit.y -= cal_shift


def main(filename, data_type):
    start_time = time.time()

    if data_type == 'data':
        file = TFile.Open('./tb16_cd_nn_reg9_nocm_corr_wfita_reco/' + filename)
        tree = file.apv_reco
        output_file = TFile('./result_trees/extracted_' + filename, 'recreate')
        output_tree = TTree('data', 'Extracted Data')
    elif data_type == 'mc':
        file = TFile.Open('./trees_raw/T16NST5G_22_03-11_16outputfile.root')
        tree = file.Lcal
        output_file = TFile('extracted_mc.root', 'recreate')
        output_tree = TTree('mc', 'Extracted MC')

    n_events = tree.GetEntries()

    tr1_n_hits = array.array('i', [0])
    tr1_hit_pad = array.array('i', [0] * 128)
    tr1_hit_sector = array.array('i', [0] * 128)
    tr1_hit_layer = array.array('i', [0] * 128)
    tr1_hit_rho = array.array('f', [0.0] * 128)
    tr1_hit_x = array.array('f', [0.0] * 128)
    tr1_hit_y = array.array('f', [0.0] * 128)
    tr1_hit_energy = array.array('f', [0.0] * 128)

    tr1_n_clusters = array.array('i', [0])
    tr1_cluster_pad = array.array('f', [0.0] * 128)
    tr1_cluster_sector = array.array('f', [0.0] * 128)
    tr1_cluster_layer = array.array('f', [0.0] * 128)
    tr1_cluster_rho = array.array('f', [0.0] * 128)
    tr1_cluster_x = array.array('f', [0.0] * 128)
    tr1_cluster_y = array.array('f', [0.0] * 128)
    tr1_cluster_energy = array.array('f', [0.0] * 128)
    tr1_cluster_n_pads = array.array('i', [0] * 128)

    tr2_n_hits = array.array('i', [0])
    tr2_hit_pad = array.array('i', [0] * 128)
    tr2_hit_sector = array.array('i', [0] * 128)
    tr2_hit_layer = array.array('i', [0] * 128)
    tr2_hit_rho = array.array('f', [0.0] * 128)
    tr2_hit_x = array.array('f', [0.0] * 128)
    tr2_hit_y = array.array('f', [0.0] * 128)
    tr2_hit_energy = array.array('f', [0.0] * 128)

    tr2_n_clusters = array.array('i', [0])
    tr2_cluster_pad = array.array('f', [0.0] * 128)
    tr2_cluster_sector = array.array('f', [0.0] * 128)
    tr2_cluster_layer = array.array('f', [0.0] * 128)
    tr2_cluster_rho = array.array('f', [0.0] * 128)
    tr2_cluster_x = array.array('f', [0.0] * 128)
    tr2_cluster_y = array.array('f', [0.0] * 128)
    tr2_cluster_energy = array.array('f', [0.0] * 128)
    tr2_cluster_n_pads = array.array('i', [0] * 128)

    cal_n_hits = array.array('i', [0])
    cal_hit_pad = array.array('i', [0] * 128 * 5)
    cal_hit_sector = array.array('i', [0] * 128 * 5)
    cal_hit_layer = array.array('i', [0] * 128 * 5)
    cal_hit_rho = array.array('f', [0.0] * 128 * 5)
    cal_hit_x = array.array('f', [0.0] * 128 * 5)
    cal_hit_y = array.array('f', [0.0] * 128 * 5)
    cal_hit_energy = array.array('f', [0.0] * 128 * 5)

    cal_n_towers = array.array('i', [0])
    cal_tower_pad = array.array('i', [0] * 128)
    cal_tower_sector = array.array('i', [0] * 128)
    cal_tower_energy = array.array('f', [0.0] * 128)
    cal_tower_cluster = array.array('i', [0] * 128)

    cal_n_clusters = array.array('i', [0])
    cal_cluster_pad = array.array('f', [0.0] * 128 * 5)
    cal_cluster_sector = array.array('f', [0.0] * 128 * 5)
    cal_cluster_layer = array.array('f', [0.0] * 128 * 5)
    cal_cluster_rho = array.array('f', [0.0] * 128 * 5)
    cal_cluster_x = array.array('f', [0.0] * 128 * 5)
    cal_cluster_y = array.array('f', [0.0] * 128 * 5)
    cal_cluster_energy = array.array('f', [0.0] * 128 * 5)
    cal_cluster_n_pads = array.array('i', [0] * 128 * 5)

    output_tree.Branch('tr1_n_hits', tr1_n_hits, 'tr1_n_hits/I')
    output_tree.Branch('tr1_hit_pad', tr1_hit_pad, 'tr1_hit_pad[tr1_n_hits]/I')
    output_tree.Branch('tr1_hit_sector', tr1_hit_sector, 'tr1_hit_sector[tr1_n_hits]/I')
    output_tree.Branch('tr1_hit_layer', tr1_hit_layer, 'tr1_hit_layer[tr1_n_hits]/I')
    output_tree.Branch('tr1_hit_rho', tr1_hit_rho, 'tr1_hit_rho[tr1_n_hits]/F')
    output_tree.Branch('tr1_hit_x', tr1_hit_x, 'tr1_hit_x[tr1_n_hits]/F')
    output_tree.Branch('tr1_hit_y', tr1_hit_y, 'tr1_hit_y[tr1_n_hits]/F')
    output_tree.Branch('tr1_hit_energy', tr1_hit_energy, 'tr1_hit_energy[tr1_n_hits]/F')

    output_tree.Branch('tr1_n_clusters', tr1_n_clusters, 'tr1_n_clusters/I')
    output_tree.Branch('tr1_cluster_pad', tr1_cluster_pad, 'tr1_cluster_pad[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_sector', tr1_cluster_sector, 'tr1_cluster_sector[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_layer', tr1_cluster_layer, 'tr1_cluster_layer[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_rho', tr1_cluster_rho, 'tr1_cluster_rho[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_x', tr1_cluster_x, 'tr1_cluster_x[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_y', tr1_cluster_y, 'tr1_cluster_y[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_energy', tr1_cluster_energy, 'tr1_cluster_energy[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_n_pads', tr1_cluster_n_pads, 'tr1_cluster_n_pads[tr1_n_clusters]/I')

    output_tree.Branch('tr2_n_hits', tr2_n_hits, 'tr2_n_hits/I')
    output_tree.Branch('tr2_hit_pad', tr2_hit_pad, 'tr2_hit_pad[tr2_n_hits]/I')
    output_tree.Branch('tr2_hit_sector', tr2_hit_sector, 'tr2_hit_sector[tr2_n_hits]/I')
    output_tree.Branch('tr2_hit_layer', tr2_hit_layer, 'tr2_hit_layer[tr2_n_hits]/I')
    output_tree.Branch('tr2_hit_rho', tr2_hit_rho, 'tr2_hit_rho[tr2_n_hits]/F')
    output_tree.Branch('tr2_hit_x', tr2_hit_x, 'tr2_hit_x[tr2_n_hits]/F')
    output_tree.Branch('tr2_hit_y', tr2_hit_y, 'tr2_hit_y[tr2_n_hits]/F')
    output_tree.Branch('tr2_hit_energy', tr2_hit_energy, 'tr2_hit_energy[tr2_n_hits]/F')

    output_tree.Branch('tr2_n_clusters', tr2_n_clusters, 'tr2_n_clusters/I')
    output_tree.Branch('tr2_cluster_pad', tr2_cluster_pad, 'tr2_cluster_pad[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_sector', tr2_cluster_sector, 'tr2_cluster_sector[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_layer', tr2_cluster_layer, 'tr2_cluster_layer[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_rho', tr2_cluster_rho, 'tr2_cluster_rho[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_x', tr2_cluster_x, 'tr2_cluster_x[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_y', tr2_cluster_y, 'tr2_cluster_y[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_energy', tr2_cluster_energy, 'tr2_cluster_energy[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_n_pads', tr2_cluster_n_pads, 'tr2_cluster_n_pads[tr2_n_clusters]/I')

    output_tree.Branch('cal_n_hits', cal_n_hits, 'cal_n_hits/I')
    output_tree.Branch('cal_hit_pad', cal_hit_pad, 'cal_hit_pad[cal_n_hits]/I')
    output_tree.Branch('cal_hit_sector', cal_hit_sector, 'cal_hit_sector[cal_n_hits]/I')
    output_tree.Branch('cal_hit_layer', cal_hit_layer, 'cal_hit_layer[cal_n_hits]/I')
    output_tree.Branch('cal_hit_rho', cal_hit_rho, 'cal_hit_rho[cal_n_hits]/F')
    output_tree.Branch('cal_hit_x', cal_hit_x, 'cal_hit_x[cal_n_hits]/F')
    output_tree.Branch('cal_hit_y', cal_hit_y, 'cal_hit_y[cal_n_hits]/F')
    output_tree.Branch('cal_hit_energy', cal_hit_energy, 'cal_hit_energy[cal_n_hits]/F')

    output_tree.Branch('cal_n_towers', cal_n_towers, 'cal_n_towers/I')
    output_tree.Branch('cal_tower_pad', cal_tower_pad, 'cal_tower_pad[cal_n_towers]/I')
    output_tree.Branch('cal_tower_sector', cal_tower_sector, 'cal_tower_sector[cal_n_towers]/I')
    output_tree.Branch('cal_tower_energy', cal_tower_energy, 'cal_tower_energy[cal_n_towers]/F')
    output_tree.Branch('cal_tower_cluster', cal_tower_cluster, 'cal_tower_cluster[cal_n_towers]/I')

    output_tree.Branch('cal_n_clusters', cal_n_clusters, 'cal_n_clusters/I')
    output_tree.Branch('cal_cluster_pad', cal_cluster_pad, 'cal_cluster_pad[cal_n_clusters]/F')
    output_tree.Branch('cal_cluster_sector', cal_cluster_sector, 'cal_cluster_sector[cal_n_clusters]/F')
    output_tree.Branch('cal_cluster_layer', cal_cluster_layer, 'cal_cluster_layer[cal_n_clusters]/F')
    output_tree.Branch('cal_cluster_rho', cal_cluster_rho, 'cal_cluster_rho[cal_n_clusters]/F')
    output_tree.Branch('cal_cluster_x', cal_cluster_x, 'cal_cluster_x[cal_n_clusters]/F')
    output_tree.Branch('cal_cluster_y', cal_cluster_y, 'cal_cluster_y[cal_n_clusters]/F')
    output_tree.Branch('cal_cluster_energy', cal_cluster_energy, 'cal_cluster_energy[cal_n_clusters]/F')
    output_tree.Branch('cal_cluster_n_pads', cal_cluster_n_pads, 'cal_cluster_n_pads[cal_n_clusters]/I')

    if data_type == 'mc':
        tr1_true_hit_rho = array.array('f', [0])
        tr1_true_hit_x = array.array('f', [0])
        tr1_true_hit_y = array.array('f', [0])
        tr2_true_hit_rho = array.array('f', [0])
        tr2_true_hit_x = array.array('f', [0])
        tr2_true_hit_y = array.array('f', [0])
        cal_true_hit_rho = array.array('f', [0])
        cal_true_hit_x = array.array('f', [0])
        cal_true_hit_y = array.array('f', [0])
        output_tree.Branch('tr1_true_hit_rho', tr1_true_hit_rho, 'tr1_true_hit_rho/F')
        output_tree.Branch('tr1_true_hit_x', tr1_true_hit_x, 'tr1_true_hit_x/F')
        output_tree.Branch('tr1_true_hit_y', tr1_true_hit_y, 'tr1_true_hit_y/F')
        output_tree.Branch('tr2_true_hit_rho', tr2_true_hit_rho, 'tr2_true_hit_rho/F')
        output_tree.Branch('tr2_true_hit_x', tr2_true_hit_x, 'tr2_true_hit_x/F')
        output_tree.Branch('tr2_true_hit_y', tr2_true_hit_y, 'tr2_true_hit_y/F')
        output_tree.Branch('cal_true_hit_rho', cal_true_hit_rho, 'cal_true_hit_rho/F')
        output_tree.Branch('cal_true_hit_x', cal_true_hit_x, 'cal_true_hit_x/F')
        output_tree.Branch('cal_true_hit_y', cal_true_hit_y, 'cal_true_hit_y/F')

    for idx, event in enumerate(tree):
        # if idx == 10:
        #    break

        if idx % (1000) == 0:
            time_min = (time.time() - start_time) // 60
            time_sec = (time.time() - start_time) % 60
            print('%d/%d events' % (idx, n_events), end=' ')
            print('%d min' % time_min, end=' ')
            print('%d sec' % time_sec)

        # Extract data.
        if data_type == 'data':
            hits_tr1, hits_tr2, hits_cal = extract_hits(event)
            align_detector(hits_tr1, hits_tr2, hits_cal)

        elif data_type == 'mc':
            hits_tr1, hits_tr2, hits_cal, true_hits = extract_mc(event)
            tr1_true_hit_x[0] = true_hits[0][0]
            tr1_true_hit_y[0] = true_hits[0][1]
            tr1_true_hit_rho[0] = true_hits[0][2]
            tr2_true_hit_x[0] = true_hits[1][0]
            tr2_true_hit_y[0] = true_hits[1][1]
            tr2_true_hit_rho[0] = true_hits[1][2]
            cal_true_hit_x[0] = true_hits[2][0]
            cal_true_hit_y[0] = true_hits[2][1]
            cal_true_hit_rho[0] = true_hits[2][2]

        towers_tr1 = set_towers(hits_tr1)
        towers_tr2 = set_towers(hits_tr2)
        towers_cal = set_towers(hits_cal)

        clusters_tr1 = set_clusters(towers_tr1, det='Tr1')
        clusters_tr2 = set_clusters(towers_tr2, det='Tr2')
        clusters_cal = set_clusters(towers_cal, det='Cal')


        if len(clusters_cal) != 0:
            clusters_tr1.sort(key=lambda x: abs(x.rho - clusters_cal[0].rho))
            clusters_tr2.sort(key=lambda x: abs(x.rho - clusters_cal[0].rho))

        tr1_n_hits[0] = len(hits_tr1)
        for i, hit in enumerate(hits_tr1):
            tr1_hit_pad[i] = hit.pad
            tr1_hit_sector[i] = hit.sector
            tr1_hit_layer[i] = hit.layer
            tr1_hit_rho[i] = hit.rho
            tr1_hit_x[i] = hit.x
            tr1_hit_y[i] = hit.y
            tr1_hit_energy[i] = hit.energy

        tr1_n_clusters[0] = len(clusters_tr1)
        for i, cluster in enumerate(clusters_tr1):
            tr1_cluster_pad[i] = cluster.pad
            tr1_cluster_sector[i] = cluster.sector
            tr1_cluster_layer[i] = cluster.layer
            tr1_cluster_rho[i] = cluster.rho
            tr1_cluster_x[i] = cluster.x
            tr1_cluster_y[i] = cluster.y
            tr1_cluster_energy[i] = cluster.energy
            tr1_cluster_n_pads[i] = cluster.n_pads

        tr2_n_hits[0] = len(hits_tr2)
        for i, hit in enumerate(hits_tr2):
            tr2_hit_pad[i] = hit.pad
            tr2_hit_sector[i] = hit.sector
            tr2_hit_layer[i] = hit.layer
            tr2_hit_rho[i] = hit.rho
            tr2_hit_x[i] = hit.x
            tr2_hit_y[i] = hit.y
            tr2_hit_energy[i] = hit.energy

        tr2_n_clusters[0] = len(clusters_tr2)
        for i, cluster in enumerate(clusters_tr2):
            tr2_cluster_pad[i] = cluster.pad
            tr2_cluster_sector[i] = cluster.sector
            tr2_cluster_layer[i] = cluster.layer
            tr2_cluster_rho[i] = cluster.rho
            tr2_cluster_x[i] = cluster.x
            tr2_cluster_y[i] = cluster.y
            tr2_cluster_energy[i] = cluster.energy
            tr2_cluster_n_pads[i] = cluster.n_pads

        cal_n_hits[0] = len(hits_cal)
        for i, hit in enumerate(hits_cal):
            cal_hit_pad[i] = hit.pad
            cal_hit_sector[i] = hit.sector
            cal_hit_layer[i] = hit.layer
            cal_hit_rho[i] = hit.rho
            cal_hit_x[i] = hit.x
            cal_hit_y[i] = hit.y
            cal_hit_energy[i] = hit.energy

        cal_n_towers[0] = len(towers_cal)
        for i, tower in enumerate(towers_cal):
            cal_tower_pad[i] = tower.pad
            cal_tower_sector[i] = tower.sector
            cal_tower_energy[i] = tower.energy
            cal_tower_cluster[i] = tower.cluster

        cal_n_clusters[0] = len(clusters_cal)
        for i, cluster in enumerate(clusters_cal):
            cal_cluster_pad[i] = cluster.pad
            cal_cluster_sector[i] = cluster.sector
            cal_cluster_layer[i] = cluster.layer
            cal_cluster_rho[i] = cluster.rho
            cal_cluster_x[i] = cluster.x
            cal_cluster_y[i] = cluster.y
            cal_cluster_energy[i] = cluster.energy
            cal_cluster_n_pads[i] = cluster.n_pads

        output_tree.Fill()

    output_tree.Write()
    output_file.Close()


main('run741_5gev.root', 'data')
main('run742_4gev.root', 'data')
main('run745_3gev.root', 'data')
main('run747_2gev.root', 'data')
main('run750_1gev.root', 'data')
# main('mc')
