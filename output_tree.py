from ROOT import TFile, TTree
import array

class OutputTree:
    def __init__(self):
        self.output_file = TFile("./extracted_mc_RENAME.root", "RECREATE")
        self.output_tree = TTree('lumical', 'MC preprocessed')
        self.define_variables()
        self.define_leafs()

    def define_variables(self):
        self.n_triggers = array.array('i', [0])
        self.trigger1 = array.array('i', [0])
        self.trigger2 = array.array('i', [0])
        self.trigger3 = array.array('i', [0])

        self.tr1_n_hits = array.array('i', [0])
        self.tr1_hit_pad = array.array('i', [0] * 128)
        self.tr1_hit_sector = array.array('i', [0] * 128)
        self.tr1_hit_layer = array.array('i', [0] * 128)
        self.tr1_hit_x = array.array('f', [0.0] * 128)
        self.tr1_hit_y = array.array('f', [0.0] * 128)
        self.tr1_hit_energy = array.array('f', [0.0] * 128)

        self.tr1_hit_type = array.array('i', [0] * 128)
        self.tr1_track_len = array.array('f', [0.0] * 128)
        self.tr1_particle_x = array.array('f', [0.0] * 128)
        self.tr1_particle_y = array.array('f', [0.0] * 128)
        self.tr1_particle_z = array.array('f', [0.0] * 128)
        self.tr1_particle_px = array.array('f', [0.0] * 128)
        self.tr1_particle_py = array.array('f', [0.0] * 128)
        self.tr1_particle_pz = array.array('f', [0.0] * 128)
        self.tr1_particle_energy = array.array('f', [0.0] * 128)

        self.tr1_n_clusters = array.array('i', [0])
        self.tr1_cluster_n_pads = array.array('i', [0] * 128)
        self.tr1_cluster_pad = array.array('f', [0.0] * 128)
        self.tr1_cluster_sector = array.array('f', [0.0] * 128)
        self.tr1_cluster_x = array.array('f', [0.0] * 128)
        self.tr1_cluster_y = array.array('f', [0.0] * 128)
        self.tr1_cluster_energy = array.array('f', [0.0] * 128)

        self.tr2_n_hits = array.array('i', [0])
        self.tr2_hit_pad = array.array('i', [0] * 128)
        self.tr2_hit_sector = array.array('i', [0] * 128)
        self.tr2_hit_layer = array.array('i', [0] * 128)
        self.tr2_hit_x = array.array('f', [0.0] * 128)
        self.tr2_hit_y = array.array('f', [0.0] * 128)
        self.tr2_hit_energy = array.array('f', [0.0] * 128)

        self.tr2_hit_type = array.array('i', [0] * 128)
        self.tr2_track_len = array.array('f', [0.0] * 128)
        self.tr2_particle_x = array.array('f', [0.0] * 128)
        self.tr2_particle_y = array.array('f', [0.0] * 128)
        self.tr2_particle_z = array.array('f', [0.0] * 128)
        self.tr2_particle_px = array.array('f', [0.0] * 128)
        self.tr2_particle_py = array.array('f', [0.0] * 128)
        self.tr2_particle_pz = array.array('f', [0.0] * 128)
        self.tr2_particle_energy = array.array('f', [0.0] * 128)

        self.tr2_n_clusters = array.array('i', [0])
        self.tr2_cluster_n_pads = array.array('i', [0] * 128)
        self.tr2_cluster_pad = array.array('f', [0.0] * 128)
        self.tr2_cluster_sector = array.array('f', [0.0] * 128)
        self.tr2_cluster_x = array.array('f', [0.0] * 128)
        self.tr2_cluster_y = array.array('f', [0.0] * 128)
        self.tr2_cluster_energy = array.array('f', [0.0] * 128)

        self.cal_n_hits = array.array('i', [0])
        self.cal_hit_pad = array.array('i', [0] * 128 * 5)
        self.cal_hit_sector = array.array('i', [0] * 128 * 5)
        self.cal_hit_layer = array.array('i', [0] * 128 * 5)
        self.cal_hit_x = array.array('f', [0.0] * 128 * 5)
        self.cal_hit_y = array.array('f', [0.0] * 128 * 5)
        self.cal_hit_energy = array.array('f', [0.0] * 128 * 5)

        self.cal_n_clusters = array.array('i', [0])
        self.cal_cluster_n_pads = array.array('i', [0] * 128 * 5)
        self.cal_cluster_n_towers = array.array('i', [0] * 128 * 5)
        self.cal_cluster_pad = array.array('f', [0.0] * 128 * 5)
        self.cal_cluster_sector = array.array('f', [0.0] * 128 * 5)
        self.cal_cluster_layer = array.array('f', [0.0] * 128 * 5)
        self.cal_cluster_x = array.array('f', [0.0] * 128 * 5)
        self.cal_cluster_y = array.array('f', [0.0] * 128 * 5)
        self.cal_cluster_energy = array.array('f', [0.0] * 128 * 5)

    def define_leafs(self):
        self.output_tree.Branch('n_triggers', self.n_triggers, 'n_triggers/I')
        self.output_tree.Branch('trigger1', self.trigger1, 'trigger1/I')
        self.output_tree.Branch('trigger2', self.trigger2, 'trigger2/I')
        self.output_tree.Branch('trigger3', self.trigger3, 'trigger3/I')

        self.output_tree.Branch('tr1_n_hits', self.tr1_n_hits, 'tr1_n_hits/I')
        self.output_tree.Branch('tr1_hit_pad', self.tr1_hit_pad, 'tr1_hit_pad[tr1_n_hits]/I')
        self.output_tree.Branch('tr1_hit_sector', self.tr1_hit_sector, 'tr1_hit_sector[tr1_n_hits]/I')
        self.output_tree.Branch('tr1_hit_layer', self.tr1_hit_layer, 'tr1_hit_layer[tr1_n_hits]/I')
        self.output_tree.Branch('tr1_hit_x', self.tr1_hit_x, 'tr1_hit_x[tr1_n_hits]/F')
        self.output_tree.Branch('tr1_hit_y', self.tr1_hit_y, 'tr1_hit_y[tr1_n_hits]/F')
        self.output_tree.Branch('tr1_hit_energy', self.tr1_hit_energy, 'tr1_hit_energy[tr1_n_hits]/F')

        self.output_tree.Branch('tr1_hit_type', self.tr1_hit_type, 'tr1_hit_type[tr1_n_hits]/I')
        self.output_tree.Branch('tr1_track_len', self.tr1_track_len, 'tr1_track_len[tr1_n_hits]/F')
        self.output_tree.Branch('tr1_particle_x', self.tr1_particle_x, 'tr1_particle_x[tr1_n_hits]/F')
        self.output_tree.Branch('tr1_particle_y', self.tr1_particle_y, 'tr1_particle_y[tr1_n_hits]/F')
        self.output_tree.Branch('tr1_particle_z', self.tr1_particle_z, 'tr1_particle_z[tr1_n_hits]/F')
        self.output_tree.Branch('tr1_particle_px', self.tr1_particle_px, 'tr1_particle_px[tr1_n_hits]/F')
        self.output_tree.Branch('tr1_particle_py', self.tr1_particle_py, 'tr1_particle_py[tr1_n_hits]/F')
        self.output_tree.Branch('tr1_particle_pz', self.tr1_particle_pz, 'tr1_particle_pz[tr1_n_hits]/F')
        self.output_tree.Branch('tr1_particle_energy', self.tr1_particle_energy, 'tr1_particle_energy[tr1_n_hits]/F')

        self.output_tree.Branch('tr1_n_clusters', self.tr1_n_clusters, 'tr1_n_clusters/I')
        self.output_tree.Branch('tr1_cluster_n_pads', self.tr1_cluster_n_pads, 'tr1_cluster_n_pads[tr1_n_clusters]/I')
        self.output_tree.Branch('tr1_cluster_pad', self.tr1_cluster_pad, 'tr1_cluster_pad[tr1_n_clusters]/F')
        self.output_tree.Branch('tr1_cluster_sector', self.tr1_cluster_sector, 'tr1_cluster_sector[tr1_n_clusters]/F')
        self.output_tree.Branch('tr1_cluster_x', self.tr1_cluster_x, 'tr1_cluster_x[tr1_n_clusters]/F')
        self.output_tree.Branch('tr1_cluster_y', self.tr1_cluster_y, 'tr1_cluster_y[tr1_n_clusters]/F')
        self.output_tree.Branch('tr1_cluster_energy', self.tr1_cluster_energy, 'tr1_cluster_energy[tr1_n_clusters]/F')

        self.output_tree.Branch('tr2_n_hits', self.tr2_n_hits, 'tr2_n_hits/I')
        self.output_tree.Branch('tr2_hit_pad', self.tr2_hit_pad, 'tr2_hit_pad[tr2_n_hits]/I')
        self.output_tree.Branch('tr2_hit_sector', self.tr2_hit_sector, 'tr2_hit_sector[tr2_n_hits]/I')
        self.output_tree.Branch('tr2_hit_layer', self.tr2_hit_layer, 'tr2_hit_layer[tr2_n_hits]/I')
        self.output_tree.Branch('tr2_hit_x', self.tr2_hit_x, 'tr2_hit_x[tr2_n_hits]/F')
        self.output_tree.Branch('tr2_hit_y', self.tr2_hit_y, 'tr2_hit_y[tr2_n_hits]/F')
        self.output_tree.Branch('tr2_hit_energy', self.tr2_hit_energy, 'tr2_hit_energy[tr2_n_hits]/F')

        self.output_tree.Branch('tr2_hit_type', self.tr2_hit_type, 'tr2_hit_type[tr2_n_hits]/I')
        self.output_tree.Branch('tr2_track_len', self.tr2_track_len, 'tr2_track_len[tr2_n_hits]/F')
        self.output_tree.Branch('tr2_particle_x', self.tr2_particle_x, 'tr2_particle_x[tr2_n_hits]/F')
        self.output_tree.Branch('tr2_particle_y', self.tr2_particle_y, 'tr2_particle_y[tr2_n_hits]/F')
        self.output_tree.Branch('tr2_particle_z', self.tr2_particle_z, 'tr2_particle_z[tr2_n_hits]/F')
        self.output_tree.Branch('tr2_particle_px', self.tr2_particle_px, 'tr2_particle_px[tr2_n_hits]/F')
        self.output_tree.Branch('tr2_particle_py', self.tr2_particle_py, 'tr2_particle_py[tr2_n_hits]/F')
        self.output_tree.Branch('tr2_particle_pz', self.tr2_particle_pz, 'tr2_particle_pz[tr2_n_hits]/F')
        self.output_tree.Branch('tr2_particle_energy', self.tr2_particle_energy, 'tr2_particle_energy[tr2_n_hits]/F')

        self.output_tree.Branch('tr2_n_clusters', self.tr2_n_clusters, 'tr2_n_clusters/I')
        self.output_tree.Branch('tr2_cluster_n_pads', self.tr2_cluster_n_pads, 'tr2_cluster_n_pads[tr2_n_clusters]/I')
        self.output_tree.Branch('tr2_cluster_pad', self.tr2_cluster_pad, 'tr2_cluster_pad[tr2_n_clusters]/F')
        self.output_tree.Branch('tr2_cluster_sector', self.tr2_cluster_sector, 'tr2_cluster_sector[tr2_n_clusters]/F')
        self.output_tree.Branch('tr2_cluster_x', self.tr2_cluster_x, 'tr2_cluster_x[tr2_n_clusters]/F')
        self.output_tree.Branch('tr2_cluster_y', self.tr2_cluster_y, 'tr2_cluster_y[tr2_n_clusters]/F')
        self.output_tree.Branch('tr2_cluster_energy', self.tr2_cluster_energy, 'tr2_cluster_energy[tr2_n_clusters]/F')

        self.output_tree.Branch('cal_n_hits', self.cal_n_hits, 'cal_n_hits/I')
        self.output_tree.Branch('cal_hit_pad', self.cal_hit_pad, 'cal_hit_pad[cal_n_hits]/I')
        self.output_tree.Branch('cal_hit_sector', self.cal_hit_sector, 'cal_hit_sector[cal_n_hits]/I')
        self.output_tree.Branch('cal_hit_layer', self.cal_hit_layer, 'cal_hit_layer[cal_n_hits]/I')
        self.output_tree.Branch('cal_hit_x', self.cal_hit_x, 'cal_hit_x[cal_n_hits]/F')
        self.output_tree.Branch('cal_hit_y', self.cal_hit_y, 'cal_hit_y[cal_n_hits]/F')
        self.output_tree.Branch('cal_hit_energy', self.cal_hit_energy, 'cal_hit_energy[cal_n_hits]/F')
        self.output_tree.Branch('cal_n_clusters', self.cal_n_clusters, 'cal_n_clusters/I')

        self.output_tree.Branch('cal_cluster_n_pads', self.cal_cluster_n_pads, 'cal_cluster_n_pads[cal_n_clusters]/I')
        self.output_tree.Branch('cal_cluster_n_towers', self.cal_cluster_n_towers, 'cal_cluster_n_towers[cal_n_clusters]/I')
        self.output_tree.Branch('cal_cluster_pad', self.cal_cluster_pad, 'cal_cluster_pad[cal_n_clusters]/F')
        self.output_tree.Branch('cal_cluster_sector', self.cal_cluster_sector, 'cal_cluster_sector[cal_n_clusters]/F')
        self.output_tree.Branch('cal_cluster_layer', self.cal_cluster_layer, 'cal_cluster_layer[cal_n_clusters]/F')
        self.output_tree.Branch('cal_cluster_x', self.cal_cluster_x, 'cal_cluster_x[cal_n_clusters]/F')
        self.output_tree.Branch('cal_cluster_y', self.cal_cluster_y, 'cal_cluster_y[cal_n_clusters]/F')
        self.output_tree.Branch('cal_cluster_energy', self.cal_cluster_energy, 'cal_cluster_energy[cal_n_clusters]/F')

    def fill_output_tree(self, event, tr1_hits, tr2_hits, cal_hits, tr1_clusters, tr2_clusters, cal_clusters):
        self.n_triggers[0] = event.n_triggers
        self.trigger1[0] = event.trigger1
        self.trigger2[0] = event.trigger2
        self.trigger3[0] = event.trigger3

        self.tr1_n_hits[0] = len(tr1_hits)
        for i, hit in enumerate(tr1_hits):
            self.tr1_hit_pad[i] = hit.pad
            self.tr1_hit_sector[i] = hit.sector
            self.tr1_hit_layer[i] = hit.layer
            self.tr1_hit_energy[i] = hit.energy
            self.tr1_hit_x[i] = hit.x
            self.tr1_hit_y[i] = hit.y

            self.tr1_hit_type[i] = hit.type
            self.tr1_track_len[i] = hit.track_len
            self.tr1_particle_x[i] = hit.p_x
            self.tr1_particle_y[i] = hit.p_y
            self.tr1_particle_z[i] = hit.p_z
            self.tr1_particle_px[i] = hit.p_px
            self.tr1_particle_py[i] = hit.p_py
            self.tr1_particle_pz[i] = hit.p_pz
            self.tr1_particle_energy[i] = hit.p_energy

        self.tr1_n_clusters[0] = len(tr1_clusters)
        for i, cluster in enumerate(tr1_clusters):
            self.tr1_cluster_n_pads[i] = cluster.n_pads
            self.tr1_cluster_pad[i] = cluster.pad
            self.tr1_cluster_sector[i] = cluster.sector
            self.tr1_cluster_x[i] = cluster.x
            self.tr1_cluster_y[i] = cluster.y
            self.tr1_cluster_energy[i] = cluster.energy

        self.tr2_n_hits[0] = len(tr2_hits)
        for i, hit in enumerate(tr2_hits):
            self.tr2_hit_pad[i] = hit.pad
            self.tr2_hit_sector[i] = hit.sector
            self.tr2_hit_layer[i] = hit.layer
            self.tr2_hit_energy[i] = hit.energy
            self.tr2_hit_x[i] = hit.x
            self.tr2_hit_y[i] = hit.y

            self.tr2_hit_type[i] = hit.type
            self.tr2_track_len[i] = hit.track_len
            self.tr2_particle_x[i] = hit.p_x
            self.tr2_particle_y[i] = hit.p_y
            self.tr2_particle_z[i] = hit.p_z
            self.tr2_particle_px[i] = hit.p_px
            self.tr2_particle_py[i] = hit.p_py
            self.tr2_particle_pz[i] = hit.p_pz
            self.tr2_particle_energy[i] = hit.p_energy

        self.tr2_n_clusters[0] = len(tr2_clusters)
        for i, cluster in enumerate(tr2_clusters):
            self.tr2_cluster_n_pads[i] = cluster.n_pads
            self.tr2_cluster_pad[i] = cluster.pad
            self.tr2_cluster_sector[i] = cluster.sector
            self.tr2_cluster_x[i] = cluster.x
            self.tr2_cluster_y[i] = cluster.y
            self.tr2_cluster_energy[i] = cluster.energy

        self.cal_n_hits[0] = len(cal_hits)
        for i, hit in enumerate(cal_hits):
            self.cal_hit_pad[i] = hit.pad
            self.cal_hit_sector[i] = hit.sector
            self.cal_hit_layer[i] = hit.layer
            self.cal_hit_x[i] = hit.x
            self.cal_hit_y[i] = hit.y
            self.cal_hit_energy[i] = hit.energy

        self.cal_n_clusters[0] = len(cal_clusters)
        for i, cluster in enumerate(cal_clusters):
            self.cal_cluster_n_pads[i] = cluster.n_pads
            self.cal_cluster_n_towers[i] = cluster.n_towers
            self.cal_cluster_pad[i] = cluster.pad
            self.cal_cluster_sector[i] = cluster.sector
            self.cal_cluster_layer[i] = cluster.layer
            self.cal_cluster_x[i] = cluster.x
            self.cal_cluster_y[i] = cluster.y
            self.cal_cluster_energy[i] = cluster.energy

        self.output_tree.Fill()

    def write_file(self):
        self.output_tree.Write()
        self.output_file.Close()
