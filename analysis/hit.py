import numpy as np

class Hit:
    def __init__(self, s, p, l, e_in_mev):
        self.sector = s
        self.pad = p
        self.layer = l
        self.energy = e_in_mev
        self.rho = 80. + 0.9 + 1.8 * p
        self.phi = np.pi * (0.5 + 1. / 12 - 1. / 48 - 1. / 24 * s)
        self.x = self.rho * np.cos(self.phi)
        self.y = self.rho * np.sin(self.phi)


class TrHit(Hit):
    def __init__(self, s, p, l, e_in_mev,
                 hit_type, track_len,
                 p_x, p_y, p_z, p_px, p_py, p_pz, p_energy):
        super().__init__(s, p, l, e_in_mev)

        self.type = hit_type
        self.track_len = track_len
        self.p_x = p_x
        self.p_y = p_y
        self.p_z = p_z
        self.p_px = p_px
        self.p_py = p_py
        self.p_pz = p_pz
        self.p_energy = p_energy
        self.seed = -1
        self.shared = False
        self.check_boundary()

    def check_boundary(self):
        # Get the true position in the detector coord system
        misalignment = 0
        p_x = -self.p_x
        p_y = self.p_y + 164.3 + misalignment
        rho = np.sqrt(p_x**2 + p_y**2)
        phi = np.arctan(p_y/p_x)
        print("{} < {} < {} - ? {}".format(self.phi - np.pi/12, phi, self.phi + np.pi/12, self.phi - np.pi/12 < phi < self.phi + np.pi/12))
        print("{} < {} < {} - ? {}".format(self.rho - 0.9, rho, self.rho + 0.9, self.rho - 0.9 < rho < self.rho + 0.9))
        input("wait")
