# File are taken from Sasha Borysov directory:
# /data/alzta/aborysov/tb_2016_data/code/tb16_reco_cd_tr_nn_wfita

from ROOT import TFile, TTree, TChain
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
           or (hit.layer >= 2 and (hit.energy < 0 or apv_nn_output[i] < 0.5))
           or bad_pad(hit.sector, hit.pad, hit.layer)):
            continue

        if hit.layer == 0:
            hits_tracker1.append(hit)
        elif hit.layer == 1:
            hits_tracker2.append(hit)
        else:
            hits_calorimeter.append(hit)

    return hits_tracker1, hits_tracker2, hits_calorimeter

def align_detector(hits_tr1, hits_tr2, hits_cal):
    # For runs > 734
    # tr1_shift = -0.16287581540422025
    # tr2_shift = 0.9707477456103106
    # cal_shift = -0.807871930206062
    tr1_shift = 0.
    tr2_shift = 0.
    cal_shift = 0.

    for hit in hits_tr1:
        hit.y -= tr1_shift
    for hit in hits_tr2:
        hit.y -= tr2_shift
    for hit in hits_cal:
        hit.y -= cal_shift


def main(beam_energy, run_type):
    start_time = time.time()

    tree = TChain("apv_reco")
    if run_type == 1:
        if beam_energy == 5:
            tree.Add("./tb16_cd_nn_reg9_nocm_corr_wfita_reco/run737_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root")
            tree.Add("./tb16_cd_nn_reg9_nocm_corr_wfita_reco/run738_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root")
            tree.Add("./tb16_cd_nn_reg9_nocm_corr_wfita_reco/run739_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root")
            tree.Add("./tb16_cd_nn_reg9_nocm_corr_wfita_reco/run740_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root")
            tree.Add("./tb16_cd_nn_reg9_nocm_corr_wfita_reco/run741_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root")
        elif beam_energy == 4:
            tree.Add("./tb16_cd_nn_reg9_nocm_corr_wfita_reco/run742_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root")
            tree.Add("./tb16_cd_nn_reg9_nocm_corr_wfita_reco/run743_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root")
            tree.Add("./tb16_cd_nn_reg9_nocm_corr_wfita_reco/run744_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root")
        elif beam_energy == 3:
            tree.Add("./tb16_cd_nn_reg9_nocm_corr_wfita_reco/run745_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root")
            tree.Add("./tb16_cd_nn_reg9_nocm_corr_wfita_reco/run746_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root")
        elif beam_energy == 2:
            tree.Add("./tb16_cd_nn_reg9_nocm_corr_wfita_reco/run747_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root")
            tree.Add("./tb16_cd_nn_reg9_nocm_corr_wfita_reco/run748_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root")
        elif beam_energy == 1:
            tree.Add("./tb16_cd_nn_reg9_nocm_corr_wfita_reco/run749_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root")
            tree.Add("./tb16_cd_nn_reg9_nocm_corr_wfita_reco/run750_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root")

    elif run_type == 2:
        if beam_energy == 5:
            tree.Add("./tb16_cd_nn_reg9_nocm_corr_wfita_reco/run754_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root")
            tree.Add("./tb16_cd_nn_reg9_nocm_corr_wfita_reco/run755_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root")
            tree.Add("./tb16_cd_nn_reg9_nocm_corr_wfita_reco/run756_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root")
        elif beam_energy == 4:
            tree.Add("./tb16_cd_nn_reg9_nocm_corr_wfita_reco/run757_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root")
            tree.Add("./tb16_cd_nn_reg9_nocm_corr_wfita_reco/run758_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root")
        elif beam_energy == 3:
            tree.Add("./tb16_cd_nn_reg9_nocm_corr_wfita_reco/run759_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root")
            tree.Add("./tb16_cd_nn_reg9_nocm_corr_wfita_reco/run760_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root")
        elif beam_energy == 2:
            tree.Add("./tb16_cd_nn_reg9_nocm_corr_wfita_reco/run761_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root")
            tree.Add("./tb16_cd_nn_reg9_nocm_corr_wfita_reco/run762_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root")
        elif beam_energy == 1:
            tree.Add("./tb16_cd_nn_reg9_nocm_corr_wfita_reco/run763_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root")
            tree.Add("./tb16_cd_nn_reg9_nocm_corr_wfita_reco/run764_tb16_charge_div_nn_reg9_nocm_corr_wfita_reco.root")


        output_file = TFile('./extracted_trees/extracted_{}_gev_run_type_{}.root'.format(beam_energy, run_type), 'recreate')
        output_tree = TTree('data', 'Extracted Data')

    for idx, event in enumerate(tree):
        # if idx == 10:
        #    break

        if idx % (10000) == 0 and idx != 0:
            time_min = (time.time() - start_time) // 60
            time_sec = (time.time() - start_time) % 60
            eta_min = (time_min * 60 + time_sec) * (tree.GetEntries() // idx - 1) // 60
            eta_sec = (time_min * 60 + time_sec) * (tree.GetEntries() // idx - 1) % 60
            print('%d min %d sec' % (time_min, time_sec), end=' ')
            print('ETA: %d min %d sec' % (eta_min, eta_sec))

        # Extract data.
        hits_tr1, hits_tr2, hits_cal = extract_hits(event)
        align_detector(hits_tr1, hits_tr2, hits_cal)

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
            # tr1_hit_pad[i] = hit.pad
            # tr1_hit_sector[i] = hit.sector
            # tr1_hit_layer[i] = hit.layer
            # tr1_hit_rho[i] = hit.rho
            # tr1_hit_x[i] = hit.x
            # tr1_hit_y[i] = hit.y
            tr1_hit_energy[i] = hit.energy

        tr1_n_clusters[0] = len(clusters_tr1)
        for i, cluster in enumerate(clusters_tr1):
            # tr1_cluster_pad[i] = cluster.pad
            # tr1_cluster_sector[i] = cluster.sector
            # tr1_cluster_layer[i] = cluster.layer
            # tr1_cluster_rho[i] = cluster.rho
            # tr1_cluster_x[i] = cluster.x
            tr1_cluster_y[i] = cluster.y
            tr1_cluster_energy[i] = cluster.energy
            # tr1_cluster_n_pads[i] = cluster.n_pads

        tr2_n_hits[0] = len(hits_tr2)
        for i, hit in enumerate(hits_tr2):
            # tr2_hit_pad[i] = hit.pad
            # tr2_hit_sector[i] = hit.sector
            # tr2_hit_layer[i] = hit.layer
            # tr2_hit_rho[i] = hit.rho
            # tr2_hit_x[i] = hit.x
            # tr2_hit_y[i] = hit.y
            tr2_hit_energy[i] = hit.energy

        tr2_n_clusters[0] = len(clusters_tr2)
        for i, cluster in enumerate(clusters_tr2):
            # tr2_cluster_pad[i] = cluster.pad
            # tr2_cluster_sector[i] = cluster.sector
            # tr2_cluster_layer[i] = cluster.layer
            # tr2_cluster_rho[i] = cluster.rho
            # tr2_cluster_x[i] = cluster.x
            tr2_cluster_y[i] = cluster.y
            tr2_cluster_energy[i] = cluster.energy
            # tr2_cluster_n_pads[i] = cluster.n_pads

        # cal_n_hits[0] = len(hits_cal)
        # for i, hit in enumerate(hits_cal):
        #     cal_hit_pad[i] = hit.pad
        #     cal_hit_sector[i] = hit.sector
        #     cal_hit_layer[i] = hit.layer
        #     cal_hit_rho[i] = hit.rho
        #     cal_hit_x[i] = hit.x
        #     cal_hit_y[i] = hit.y
        #     cal_hit_energy[i] = hit.energy

        # cal_n_towers[0] = len(towers_cal)
        # for i, tower in enumerate(towers_cal):
        #     cal_tower_pad[i] = tower.pad
        #     cal_tower_sector[i] = tower.sector
        #     cal_tower_energy[i] = tower.energy
        #     cal_tower_cluster[i] = tower.cluster

        cal_n_clusters[0] = len(clusters_cal)
        for i, cluster in enumerate(clusters_cal):
        #     cal_cluster_pad[i] = cluster.pad
        #     cal_cluster_sector[i] = cluster.sector
        #     cal_cluster_layer[i] = cluster.layer
        #     cal_cluster_rho[i] = cluster.rho
        #     cal_cluster_x[i] = cluster.x
            cal_cluster_y[i] = cluster.y
        #     cal_cluster_energy[i] = cluster.energy
        #     cal_cluster_n_pads[i] = cluster.n_pads

        output_tree.Fill()

    output_tree.Write()
    output_file.Close()

    print("Hooray, extracted tree file is ready, take it :3")


# main('data', beam_energy=1, run_type=1)
# main('data', beam_energy=2, run_type=1)
# main('data', beam_energy=3, run_type=1)
# main('data', beam_energy=4, run_type=1)
# main('data', beam_energy=5, run_type=1)
main('data', beam_energy=5, run_type=3)

# main('mc')
