import numpy as np
from ROOT import TTree, TFile
import array
from scipy import special

f=open("./noise.txt", "r")
lines = f.readlines()
noise = np.zeros((4, 64, 8)) # In MeV
for line in lines[1:]:
    values = line.split('  ')
    sector = int(values[0])
    if sector == -1:
        continue
    pad = int(values[1])
    layer = int(values[2])
    noise[sector, pad, layer] = float(values[3])

f_input = TFile.Open("../trees_5gev_e/lucas_geocuts.root", 'read')
t_input = f_input.lumical

f_output = TFile.Open("../trees_5gev_e/lucas_geocuts_noise_cal_eff.root", "RECREATE")
t_output = TTree('lumical', 'MC smeared')

cal_n_hits_new = array.array('i', [0])
cal_energy_new = array.array('f', [0.0] * 128 * 5)

t_output.Branch('cal_n_hits_new', cal_n_hits_new, 'cal_n_hits_new/I')
t_output.Branch('cal_energy_new', cal_energy_new, 'cal_energy_new[cal_n_hits_new]/F')


for idx, event in enumerate(t_input):
    if (idx % 1000) == 0:
        print("Event: ", idx)
    if idx == 100000:
        break

    n_hits = event.cal_n_hits
    j = 0
    for i in range(n_hits):
        sector = event.cal_hit_sector[i]
        pad = event.cal_hit_pad[i]
        layer = event.cal_hit_layer[i]
        energy = event.cal_hit_energy[i]
        rand_noise = np.random.normal(0., 0.7*noise[sector][pad][layer])
        energy_smeared = energy + rand_noise
        if energy_smeared <= 0:
            continue
        # Itamar values
        # S0_cal = 0.819
        # p1_cal = 2.166
        S0_cal = 2.
        p1_cal = 3.3
        p0 = 0.999 / 2.
        # If hit is in calorimeter - reject hit with a probability based on efficiency curve
        # measured for the paper
        if np.random.random() > (1. + special.erf((energy_smeared/0.0885 - S0_cal) / p1_cal)) * p0:
            continue
        cal_energy_new[j] = energy_smeared
        j +=1
    cal_n_hits_new[0] = j

    t_output.Fill()

f_output.Write()
f_output.Close()
