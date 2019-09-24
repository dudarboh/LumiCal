from ROOT import TFile, TTree
import array
import numpy as np
import time
import cProfile
import pstats

# Important notes which confused me in the past
# 4 sectors: 0, 1, 2, 3
# 64 pads: 0, 1, 2, ..., 63
# 8 layers: 0, 1 - trackers; 2, 3, 4, 5, 6, 7 - calorimeter; 7 - tab (bad)


def bad_pad(sector, pad, layer):
    return ((layer == 0 and sector == 1 and pad in (26, 61))
            or (layer == 0 and sector == 2 and pad in (31, 57, 61))
            or (layer == 2 and sector == 1 and pad in (28, 31, 34))
            or (layer == 2 and sector == 2 and pad in (38, 53))
            or (layer == 3 and sector == 2 and pad in (31, 33, 52, 55, 61))
            or (layer == 4 and sector == 1 and pad in (29, 39, 41, 55, 56))
            or (layer == 4 and sector == 2 and pad in (28,))
            or (layer == 5 and sector == 1 and pad in (32, 36, 40, 41, 44, 45, 49, 56, 58))
            or (layer == 5 and sector == 2 and pad in (28, 52, 54, 61))
            or (layer == 6 and sector == 1 and pad in (26, 30))
            or (layer == 6 and sector == 2 and pad in (34, 42, 54, 57, 59, 60)))


def align_detector(hits_tr1, hits_tr2, hits_cal):
    # Was checked for Itamars misalignment.
    tr1_shift = -0.2480096316087952
    tr2_shift = 0.9012092661162399
    cal_shift = -0.6531996345075299

    for hit in hits_tr1:
        hit.y -= tr1_shift
    for hit in hits_tr2:
        hit.y -= tr2_shift
    for hit in hits_cal:
        hit.y -= cal_shift


class Hit:
    def __init__(self, sector, pad, layer, n_bs_particles, n_dir_particles, energy_in_mip):
        self.sector = sector
        self.pad = pad
        self.layer = layer
        self.energy = energy_in_mip
        self.n_bs_particles = n_bs_particles
        self.n_dir_particles = n_dir_particles
        rho = 80. + 0.9 + 1.8 * self.pad
        phi = np.pi / 2 + np.pi / 12 - np.pi / 48 - np.pi / 24 * self.sector
        self.x = rho * np.cos(phi)
        self.y = rho * np.sin(phi)


class Tower:
    def __init__(self, tower_hits):
        self.hits = tower_hits

        self.sector = self.hits[0].sector
        self.pad = self.hits[0].pad

        self.energy = sum(hit.energy for hit in self.hits)

        self.neighbor = -1
        self.cluster = -1


class Cluster:
    def __init__(self, cluster_hits, det):
        self.hits = cluster_hits
        self.det = det

        self.energy = self.get_energy()

        # Energies in trackers. LogW in Calorimeter
        self.weights = self.get_weights()
        self.sum_of_weights = sum(self.weights)
        if self.sum_of_weights != 0:
            self.sector = self.get_cluster_pos("sector")
            self.pad = self.get_cluster_pos("pad")
            self.layer = self.get_cluster_pos("layer")
            self.x = self.get_cluster_pos("x")
            self.y = self.get_cluster_pos("y")
        else:
            self.sector = -999
            self.pad = -999
            self.layer = -999
            self.x = -999
            self.y = -999

    def get_energy(self):
        """Return total energy of all hits"""
        return sum(hit.energy for hit in self.hits)

    def get_weights(self):
        """Return list of hits' weights"""
        if self.det == 'Cal':
            w0 = 3.4
            return [max(0, w0 + np.log(hit.energy / self.energy)) for hit in self.hits]
        else:
            return [hit.energy for hit in self.hits]

    def get_cluster_pos(self, pos):
        """Return weighted cluster position coordinate: x, y, pad, sector, layer"""
        return sum(getattr(hit, pos) * weight for hit, weight in zip(self.hits, self.weights)) / self.sum_of_weights

    def merge(self, cluster2):
        self.hits.extend(cluster2.hits)
        self.energy = self.get_energy()
        self.weights = self.get_weights()
        self.sum_of_weights = sum(self.weights)
        if self.sum_of_weights != 0:
            self.sector = self.get_cluster_pos("sector")
            self.pad = self.get_cluster_pos("pad")
            self.layer = self.get_cluster_pos("layer")
            self.x = self.get_cluster_pos("x")
            self.y = self.get_cluster_pos("y")
        else:
            self.sector = -999
            self.pad = -999
            self.layer = -999
            self.x = -999
            self.y = -999


def set_most_energetic_neighbors(towers_list):
    """Assign most energetic neighbor to the towers"""
    for tower in towers_list:
        neighbor = tower
        for tower_neighbor in towers_list:
            if (tower_neighbor.sector in range(tower.sector - 1, tower.sector + 2)
               and tower_neighbor.pad in range(tower.pad - 1, tower.pad + 2)
               and tower_neighbor.energy >= neighbor.energy):
                neighbor = tower_neighbor
        tower.neighbor = neighbor


def set_tower_clusters(towers_list):
    """Assign cluster indices to the towers"""
    cluster_idx = 0
    for tower in towers_list:
        if (tower.neighbor.sector, tower.neighbor.pad) == (tower.sector, tower.pad):
            tower.cluster = cluster_idx
            cluster_idx += 1

    while any(tower.cluster == -1 for tower in towers_list):
        for tower in towers_list:
            if tower.cluster == -1 and tower.neighbor.cluster != -1:
                tower.cluster = tower.neighbor.cluster


def merge_clusters(clusters_list):
    """Merge pair of clusters if they meet following condition"""
    while True:
        for clst1, clst2 in zip(clusters_list, clusters_list):
            if clst1 == clst2:
                continue
            distance = abs(clst1.y - clst2.y)
            ratio = clst2.energy / clst1.energy
            if distance < 9. or (distance < 36. and ratio < 0.1 - 0.1 / 20 * distance):
                clst1.merge(clst2)
                clusters_list.remove(clst2)
                break
        else:
            break


def make_hits_lists(event):
    hits_tracker1 = []
    hits_tracker2 = []
    hits_calorimeter = []

    for sector, pad, layer, energy, n_bs, n_dir in zip(event.hit_sector,
                                                       event.hit_pad,
                                                       event.hit_layer,
                                                       event.hit_energy,
                                                       event.n_bs_particles,
                                                       event.n_dir_particles):

        # Selection as in data
        if ((layer > 1 and energy < 1.4)
            # or pad < 20
            or layer == 7
            # or sector == 0 or sector == 3
            or (layer <= 1 and energy <= 0.)
            or bad_pad(sector, pad, layer)):
            continue

        # If passed the selection create hit and add to corresponding list
        hit = Hit(sector, pad, layer, n_bs, n_dir, energy)

        if layer == 0:
            hits_tracker1.append(hit)
        elif layer == 1:
            hits_tracker2.append(hit)
        else:
            hits_calorimeter.append(hit)

    return hits_tracker1, hits_tracker2, hits_calorimeter


def make_towers_list(hits_list):
    """Return tower list out of hits list"""
    towers_pos = set((hit.sector, hit.pad) for hit in hits_list)
    towers = []
    for pos in towers_pos:
        tower_hits = [hit for hit in hits_list if (hit.sector, hit.pad) == pos]
        towers.append(Tower(tower_hits))

    towers.sort(key=lambda x: x.energy, reverse=True)

    return towers


def make_clusters_list(towers_list, det):
    """Return cluster list out of towers list"""
    if not towers_list:
        return []
    clusters = []

    n_clusters = max(tower.cluster for tower in towers_list) + 1
    for i in range(n_clusters):
        cluster_hits = []
        for tower in towers_list:
            if tower.cluster == i:
                cluster_hits.extend(tower.hits)
        clusters.append(Cluster(cluster_hits, det))

    # Sort to start merging the most energetic ones
    clusters.sort(key=lambda x: x.energy, reverse=True)

    merge_clusters(clusters)
    clusters.sort(key=lambda x: x.energy, reverse=True)

    return clusters


def main(filename):
    start_time = time.time()
    input_file = TFile.Open('../mc/' + filename)
    input_tree = input_file.LumiCal
    print("Total n events in loaded files: ", input_tree.GetEntries())

    output_file = TFile("../extracted/extracted_" + filename, "RECREATE")
    output_tree = TTree('mc', 'Extracted MC')

    tr1_n_hits = array.array('i', [0])
    tr1_hit_pad = array.array('i', [0] * 128)
    tr1_hit_sector = array.array('i', [0] * 128)
    tr1_hit_layer = array.array('i', [0] * 128)
    tr1_hit_x = array.array('f', [0.0] * 128)
    tr1_hit_y = array.array('f', [0.0] * 128)
    tr1_hit_n_bs_particles = array.array('i', [0] * 128)
    tr1_hit_n_dir_particles = array.array('i', [0] * 128)
    tr1_hit_energy = array.array('f', [0.0] * 128)
    tr1_n_clusters = array.array('i', [0])
    tr1_cluster_pad = array.array('f', [0.0] * 128)
    tr1_cluster_sector = array.array('f', [0.0] * 128)
    tr1_cluster_layer = array.array('f', [0.0] * 128)
    tr1_cluster_x = array.array('f', [0.0] * 128)
    tr1_cluster_y = array.array('f', [0.0] * 128)
    tr1_cluster_energy = array.array('f', [0.0] * 128)

    tr2_n_hits = array.array('i', [0])
    tr2_hit_pad = array.array('i', [0] * 128)
    tr2_hit_sector = array.array('i', [0] * 128)
    tr2_hit_layer = array.array('i', [0] * 128)
    tr2_hit_x = array.array('f', [0.0] * 128)
    tr2_hit_y = array.array('f', [0.0] * 128)
    tr2_hit_n_bs_particles = array.array('i', [0] * 128)
    tr2_hit_n_dir_particles = array.array('i', [0] * 128)
    tr2_hit_energy = array.array('f', [0.0] * 128)
    tr2_n_clusters = array.array('i', [0])
    tr2_cluster_pad = array.array('f', [0.0] * 128)
    tr2_cluster_sector = array.array('f', [0.0] * 128)
    tr2_cluster_layer = array.array('f', [0.0] * 128)
    tr2_cluster_x = array.array('f', [0.0] * 128)
    tr2_cluster_y = array.array('f', [0.0] * 128)
    tr2_cluster_energy = array.array('f', [0.0] * 128)

    cal_n_hits = array.array('i', [0])
    cal_hit_pad = array.array('i', [0] * 128 * 5)
    cal_hit_sector = array.array('i', [0] * 128 * 5)
    cal_hit_layer = array.array('i', [0] * 128 * 5)
    cal_hit_x = array.array('f', [0.0] * 128 * 5)
    cal_hit_y = array.array('f', [0.0] * 128 * 5)
    cal_hit_n_bs_particles = array.array('i', [0] * 128)
    cal_hit_n_dir_particles = array.array('i', [0] * 128)
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
    cal_cluster_x = array.array('f', [0.0] * 128 * 5)
    cal_cluster_y = array.array('f', [0.0] * 128 * 5)
    cal_cluster_energy = array.array('f', [0.0] * 128 * 5)

    output_tree.Branch('tr1_n_hits', tr1_n_hits, 'tr1_n_hits/I')
    output_tree.Branch('tr1_hit_pad', tr1_hit_pad, 'tr1_hit_pad[tr1_n_hits]/I')
    output_tree.Branch('tr1_hit_sector', tr1_hit_sector, 'tr1_hit_sector[tr1_n_hits]/I')
    output_tree.Branch('tr1_hit_layer', tr1_hit_layer, 'tr1_hit_layer[tr1_n_hits]/I')
    output_tree.Branch('tr1_hit_x', tr1_hit_x, 'tr1_hit_x[tr1_n_hits]/F')
    output_tree.Branch('tr1_hit_y', tr1_hit_y, 'tr1_hit_y[tr1_n_hits]/F')
    output_tree.Branch('tr1_hit_n_bs_particles', tr1_hit_n_bs_particles, 'tr1_hit_n_bs_particles[tr1_n_hits]/I')
    output_tree.Branch('tr1_hit_n_dir_particles', tr1_hit_n_dir_particles, 'tr1_hit_n_dir_particles[tr1_n_hits]/I')
    output_tree.Branch('tr1_hit_energy', tr1_hit_energy, 'tr1_hit_energy[tr1_n_hits]/F')
    output_tree.Branch('tr1_n_clusters', tr1_n_clusters, 'tr1_n_clusters/I')
    output_tree.Branch('tr1_cluster_pad', tr1_cluster_pad, 'tr1_cluster_pad[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_sector', tr1_cluster_sector, 'tr1_cluster_sector[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_layer', tr1_cluster_layer, 'tr1_cluster_layer[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_x', tr1_cluster_x, 'tr1_cluster_x[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_y', tr1_cluster_y, 'tr1_cluster_y[tr1_n_clusters]/F')
    output_tree.Branch('tr1_cluster_energy', tr1_cluster_energy, 'tr1_cluster_energy[tr1_n_clusters]/F')

    output_tree.Branch('tr2_n_hits', tr2_n_hits, 'tr2_n_hits/I')
    output_tree.Branch('tr2_hit_pad', tr2_hit_pad, 'tr2_hit_pad[tr2_n_hits]/I')
    output_tree.Branch('tr2_hit_sector', tr2_hit_sector, 'tr2_hit_sector[tr2_n_hits]/I')
    output_tree.Branch('tr2_hit_layer', tr2_hit_layer, 'tr2_hit_layer[tr2_n_hits]/I')
    output_tree.Branch('tr2_hit_x', tr2_hit_x, 'tr2_hit_x[tr2_n_hits]/F')
    output_tree.Branch('tr2_hit_y', tr2_hit_y, 'tr2_hit_y[tr2_n_hits]/F')
    output_tree.Branch('tr2_hit_n_bs_particles', tr2_hit_n_bs_particles, 'tr2_hit_n_bs_particles[tr2_n_hits]/I')
    output_tree.Branch('tr2_hit_n_dir_particles', tr2_hit_n_dir_particles, 'tr2_hit_n_dir_particles[tr2_n_hits]/I')
    output_tree.Branch('tr2_hit_energy', tr2_hit_energy, 'tr2_hit_energy[tr2_n_hits]/F')
    output_tree.Branch('tr2_n_clusters', tr2_n_clusters, 'tr2_n_clusters/I')
    output_tree.Branch('tr2_cluster_pad', tr2_cluster_pad, 'tr2_cluster_pad[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_sector', tr2_cluster_sector, 'tr2_cluster_sector[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_layer', tr2_cluster_layer, 'tr2_cluster_layer[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_x', tr2_cluster_x, 'tr2_cluster_x[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_y', tr2_cluster_y, 'tr2_cluster_y[tr2_n_clusters]/F')
    output_tree.Branch('tr2_cluster_energy', tr2_cluster_energy, 'tr2_cluster_energy[tr2_n_clusters]/F')

    output_tree.Branch('cal_n_hits', cal_n_hits, 'cal_n_hits/I')
    output_tree.Branch('cal_hit_pad', cal_hit_pad, 'cal_hit_pad[cal_n_hits]/I')
    output_tree.Branch('cal_hit_sector', cal_hit_sector, 'cal_hit_sector[cal_n_hits]/I')
    output_tree.Branch('cal_hit_layer', cal_hit_layer, 'cal_hit_layer[cal_n_hits]/I')
    output_tree.Branch('cal_hit_x', cal_hit_x, 'cal_hit_x[cal_n_hits]/F')
    output_tree.Branch('cal_hit_y', cal_hit_y, 'cal_hit_y[cal_n_hits]/F')
    output_tree.Branch('cal_hit_n_bs_particles', cal_hit_n_bs_particles, 'cal_hit_n_bs_particles[cal_n_hits]/I')
    output_tree.Branch('cal_hit_n_dir_particles', cal_hit_n_dir_particles, 'cal_hit_n_dir_particles[cal_n_hits]/I')
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
    output_tree.Branch('cal_cluster_x', cal_cluster_x, 'cal_cluster_x[cal_n_clusters]/F')
    output_tree.Branch('cal_cluster_y', cal_cluster_y, 'cal_cluster_y[cal_n_clusters]/F')
    output_tree.Branch('cal_cluster_energy', cal_cluster_energy, 'cal_cluster_energy[cal_n_clusters]/F')

    n_events = input_tree.GetEntries()
    for idx, event in enumerate(input_tree):
        # if idx == 200000:
        #     break

        if idx % (1000) == 0:
            time_min = (time.time() - start_time) // 60
            time_sec = (time.time() - start_time) % 60
            print('Event: {} out of {};'.format(idx, n_events), end=' ')
            print('{} min {} sec'.format(time_min, time_sec))

        hits_tr1, hits_tr2, hits_cal = make_hits_lists(event)
        align_detector(hits_tr1, hits_tr2, hits_cal)

        towers_tr1 = make_towers_list(hits_tr1)
        towers_tr2 = make_towers_list(hits_tr2)
        towers_cal = make_towers_list(hits_cal)

        if towers_tr1:
            set_most_energetic_neighbors(towers_tr1)
            set_tower_clusters(towers_tr1)
        if towers_tr2:
            set_most_energetic_neighbors(towers_tr2)
            set_tower_clusters(towers_tr2)
        if towers_cal:
            set_most_energetic_neighbors(towers_cal)
            set_tower_clusters(towers_cal)

        clusters_tr1 = make_clusters_list(towers_tr1, det='Tr1')
        clusters_tr2 = make_clusters_list(towers_tr2, det='Tr2')
        clusters_cal = make_clusters_list(towers_cal, det='Cal')

        # Resort clusters in trackers by distance to main cluster in calorimeter
        if len(clusters_cal) != 0:
            clusters_tr1.sort(key=lambda x: abs(x.y - clusters_cal[0].y))
            clusters_tr2.sort(key=lambda x: abs(x.y - clusters_cal[0].y))

        tr1_n_hits[0] = len(hits_tr1)
        for i, hit in enumerate(hits_tr1):
            tr1_hit_pad[i] = hit.pad
            tr1_hit_sector[i] = hit.sector
            tr1_hit_layer[i] = hit.layer
            tr1_hit_x[i] = hit.x
            tr1_hit_y[i] = hit.y
            tr1_hit_n_bs_particles[i] = hit.n_bs_particles
            tr1_hit_n_dir_particles[i] = hit.n_dir_particles
            tr1_hit_energy[i] = hit.energy
        tr1_n_clusters[0] = len(clusters_tr1)
        for i, cluster in enumerate(clusters_tr1):
            tr1_cluster_pad[i] = cluster.pad
            tr1_cluster_sector[i] = cluster.sector
            tr1_cluster_layer[i] = cluster.layer
            tr1_cluster_x[i] = cluster.x
            tr1_cluster_y[i] = cluster.y
            tr1_cluster_energy[i] = cluster.energy

        tr2_n_hits[0] = len(hits_tr2)
        for i, hit in enumerate(hits_tr2):
            tr2_hit_pad[i] = hit.pad
            tr2_hit_sector[i] = hit.sector
            tr2_hit_layer[i] = hit.layer
            tr2_hit_x[i] = hit.x
            tr2_hit_y[i] = hit.y
            tr2_hit_n_bs_particles[i] = hit.n_bs_particles
            tr2_hit_n_dir_particles[i] = hit.n_dir_particles
            tr2_hit_energy[i] = hit.energy
        tr2_n_clusters[0] = len(clusters_tr2)
        for i, cluster in enumerate(clusters_tr2):
            tr2_cluster_pad[i] = cluster.pad
            tr2_cluster_sector[i] = cluster.sector
            tr2_cluster_layer[i] = cluster.layer
            tr2_cluster_x[i] = cluster.x
            tr2_cluster_y[i] = cluster.y
            tr2_cluster_energy[i] = cluster.energy

        cal_n_hits[0] = len(hits_cal)
        for i, hit in enumerate(hits_cal):
            cal_hit_pad[i] = hit.pad
            cal_hit_sector[i] = hit.sector
            cal_hit_layer[i] = hit.layer
            cal_hit_x[i] = hit.x
            cal_hit_y[i] = hit.y
            cal_hit_n_bs_particles[i] = hit.n_bs_particles
            cal_hit_n_dir_particles[i] = hit.n_dir_particles
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
            cal_cluster_x[i] = cluster.x
            cal_cluster_y[i] = cluster.y
            cal_cluster_energy[i] = cluster.energy

        output_tree.Fill()

    output_tree.Write()
    output_file.Close()

    print("Hooray, extracted tree file is ready, take it :3")


pr = cProfile.Profile()
pr.enable()
main('lucas_5gev.root')
main('lucas_4gev.root')
main('lucas_3gev.root')
main('lucas_2gev.root')
main('lucas_1gev.root')

pr.disable()

ps = pstats.Stats(pr).sort_stats(pstats.SortKey.CUMULATIVE)
ps.print_stats()
