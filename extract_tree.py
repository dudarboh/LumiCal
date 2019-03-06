from ROOT import TFile, TTree
import array
import time

from objects import Hit, HitMC, Tower, clustering


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

        if (hit.sector == 0 or hit.sector == 3 or hit.layer == 7
           or (hit.sector == 1 and hit.pad < 20)
           or (hit.sector == 2 and hit.pad < 20)
           or hit.sector < 0  # This one is changed due to python C++ difference in %.
           or (hit.layer < 2 and (signal_arr[i] < 0. or apv_nn_output[i] < 0.5))
           or (hit.layer >= 2 and (hit.energy < 1.4 or apv_nn_output[i] < 0.5))):
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
    z_tr1 = 3300.5109092863354
    z_tr2 = 3325.513351442981
    z_cal = 3384.032740480225

    # +177.2 is to convert y from montecarto to my coordinate system
    x_tr1 = vx + px * z_tr1 / pz
    y_tr1 = vy + py * z_tr1 / pz + 176.3 + 0.9
    rho_tr1 = (x_tr1**2 + y_tr1**2)**0.5

    x_tr2 = vx + px * z_tr2 / pz
    y_tr2 = vy + py * z_tr2 / pz + 176.3 + 0.9
    rho_tr2 = (x_tr2**2 + y_tr2**2)**0.5

    x_cal = vx + px * z_cal / pz
    y_cal = vy + py * z_cal / pz + 176.3 + 0.9
    rho_cal = (x_cal**2 + y_cal**2)**0.5

    true_hits = [(x_tr1, y_tr1, rho_tr1), (x_tr2, y_tr2, rho_tr2), (x_cal, y_cal, rho_cal)]

    n_hits = event.numHits
    hits_calorimeter = []
    hits_tracker1 = []
    hits_tracker2 = []

    for i in range(n_hits):
        cell_id = event.GetLeaf("Hits.cellID").GetValue(i)
        energy = event.GetLeaf("Hits.eHit").GetValue(i)

        hit = HitMC(cell_id, energy)

        if (hit.sector == 0 or hit.sector == 3 or layer == 7
           or (hit.sector == 1 and hit.pad < 20)
           or (hit.sector == 2 and hit.pad < 20)
           or (hit.layer >= 2 and hit.energy < 1.4)):
            continue

        if hit.layer == 0:
            hits_tracker1.append(hit)
        elif hit.layer == 1:
            hits_tracker2.append(hit)
        else:
            hits_calorimeter.append(hit)

    return hits_tracker1, hits_tracker2, hits_calorimeter, true_hits


def extract_towers(hits):
    towers_pos = set([(hit.sector, hit.pad) for hit in hits])
    towers = []
    for pos in towers_pos:
        tower_hits = [hit for hit in hits if (hit.sector, hit.pad) == pos]
        towers.append(Tower(tower_hits, 'data'))
    return towers


def main(data_type):
    start_time = time.time()

    if data_type == 'data':
        file = TFile.Open('./trees/run741_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root')
        tree = file.apv_reco
        output_file = TFile('extracted_data_RENAME.root', 'recreate')
        output_tree = TTree('data', 'Extracted Data')
    elif data_type == 'mc':
        file = TFile.Open('./trees/mc/T16NST5G_22_03-11_16outputfile.root')
        tree = file.Lcal
        output_file = TFile('extracted_mc_RENAME.root', 'recreate')
        output_tree = TTree('mc', 'Extracted MC')

    n_events = tree.GetEntries()

    tr1_n_clusters = array.array('i', [0])
    tr1_cluster_pad = array.array('f', [0.0] * 256)
    tr1_cluster_sector = array.array('f', [0.0] * 256)
    tr1_cluster_layer = array.array('f', [0.0] * 256)
    tr1_cluster_rho = array.array('f', [0.0] * 256)
    tr1_cluster_x = array.array('f', [0.0] * 256)
    tr1_cluster_y = array.array('f', [0.0] * 256)
    tr1_cluster_energy = array.array('f', [0.0] * 256)
    tr1_cluster_n_pads = array.array('f', [0.0] * 256)

    tr2_n_clusters = array.array('i', [0])
    tr2_cluster_pad = array.array('f', [0.0] * 256)
    tr2_cluster_sector = array.array('f', [0.0] * 256)
    tr2_cluster_layer = array.array('f', [0.0] * 256)
    tr2_cluster_rho = array.array('f', [0.0] * 256)
    tr2_cluster_x = array.array('f', [0.0] * 256)
    tr2_cluster_y = array.array('f', [0.0] * 256)
    tr2_cluster_energy = array.array('f', [0.0] * 256)
    tr2_cluster_n_pads = array.array('f', [0.0] * 256)

    cal_n_clusters = array.array('i', [0])
    cal_cluster_pad = array.array('f', [0.0] * 256)
    cal_cluster_sector = array.array('f', [0.0] * 256)
    cal_cluster_layer = array.array('f', [0.0] * 256)
    cal_cluster_rho = array.array('f', [0.0] * 256)
    cal_cluster_x = array.array('f', [0.0] * 256)
    cal_cluster_y = array.array('f', [0.0] * 256)
    cal_cluster_energy = array.array('f', [0.0] * 256)
    cal_cluster_n_pads = array.array('f', [0.0] * 256)

    output_tree.Branch('tr1_n_clusters', tr1_n_clusters, 'tr1_n_clusters/I')
    output_tree.Branch('tr1_cluster_pad', tr1_cluster_pad, 'tr1_cluster_pad[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_sector', tr1_cluster_sector, 'tr1_cluster_sector[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_layer', tr1_cluster_layer, 'tr1_cluster_layer[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_rho', tr1_cluster_rho, 'tr1_cluster_rho[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_x', tr1_cluster_x, 'tr1_cluster_x[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_y', tr1_cluster_y, 'tr1_cluster_y[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_energy', tr1_cluster_energy, 'tr1_cluster_energy[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_n_pads', tr1_cluster_n_pads, 'tr1_cluster_n_pads[tr1_n_clusters]/F')

    output_tree.Branch('tr2_n_clusters', tr2_n_clusters, 'tr2_n_clusters/I')
    output_tree.Branch('tr2_cluster_pad', tr2_cluster_pad, 'tr2_cluster_pad[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_sector', tr2_cluster_sector, 'tr2_cluster_sector[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_layer', tr2_cluster_layer, 'tr2_cluster_layer[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_rho', tr2_cluster_rho, 'tr2_cluster_rho[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_x', tr2_cluster_x, 'tr2_cluster_x[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_y', tr2_cluster_y, 'tr2_cluster_y[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_energy', tr2_cluster_energy, 'tr2_cluster_energy[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_n_pads', tr2_cluster_n_pads, 'tr2_cluster_n_pads[tr2_n_clusters]/F')

    output_tree.Branch('cal_n_clusters', cal_n_clusters, 'cal_n_clusters/I')
    output_tree.Branch('cal_cluster_pad', cal_cluster_pad, 'cal_cluster_pad[cal_n_clusters]/F')
    output_tree.Branch('cal_cluster_sector', cal_cluster_sector, 'cal_cluster_sector[cal_n_clusters]/F')
    output_tree.Branch('cal_cluster_layer', cal_cluster_layer, 'cal_cluster_layer[cal_n_clusters]/F')
    output_tree.Branch('cal_cluster_rho', cal_cluster_rho, 'cal_cluster_rho[cal_n_clusters]/F')
    output_tree.Branch('cal_cluster_x', cal_cluster_x, 'cal_cluster_x[cal_n_clusters]/F')
    output_tree.Branch('cal_cluster_y', cal_cluster_y, 'cal_cluster_y[cal_n_clusters]/F')
    output_tree.Branch('cal_cluster_energy', cal_cluster_energy, 'cal_cluster_energy[cal_n_clusters]/F')
    output_tree.Branch('cal_cluster_n_pads', cal_cluster_n_pads, 'cal_cluster_n_pads[cal_n_clusters]/F')

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
        # if idx != 1999:
        #     continue

        if idx % (1000) == 0:
            time_min = (time.time() - start_time) // 60
            time_sec = (time.time() - start_time) % 60
            print('%d/%d events' % (idx, n_events), end=' ')
            print('%d min' % time_min, end=' ')
            print('%d sec' % time_sec)

        # Extract data.
        if data_type == 'data':
            hits_tr1, hits_tr2, hits_cal = extract_hits(event)
        elif data_type == 'mc':
            signals_tr1, signals_tr2, signals_cal, true_hits = extract_mc(event)
            tr1_true_hit_x[0] = true_hits[0][0]
            tr1_true_hit_y[0] = true_hits[0][1]
            tr1_true_hit_rho[0] = true_hits[0][2]
            tr2_true_hit_x[0] = true_hits[1][0]
            tr2_true_hit_y[0] = true_hits[1][1]
            tr2_true_hit_rho[0] = true_hits[1][2]
            cal_true_hit_x[0] = true_hits[2][0]
            cal_true_hit_y[0] = true_hits[2][1]
            cal_true_hit_rho[0] = true_hits[2][2]

        towers_tr1 = extract_towers(hits_tr1)
        towers_tr2 = extract_towers(hits_tr2)
        towers_cal = extract_towers(hits_cal)

        clusters_tr1 = clustering(towers_tr1, merge='off', det='Tr1')
        clusters_tr2 = clustering(towers_tr2, merge='off', det='Tr2')
        clusters_cal = clustering(towers_cal, merge='on', det='Cal')

        if len(clusters_cal) != 0:
            clusters_tr1 = sorted(clusters_tr1, key=lambda x: abs(x.rho - clusters_cal[0].rho))
            clusters_tr2 = sorted(clusters_tr2, key=lambda x: abs(x.rho - clusters_cal[0].rho))

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


main('data')
