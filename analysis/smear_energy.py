import numpy as np
from ROOT import TTree, TFile
import array

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

f_output = TFile.Open("../trees_5gev_e/lucas_geocuts_noise_cal.root", "RECREATE")
t_output = TTree('lumical', 'MC smeared')

# tr1_n_hits = array.array('i', [0])
# tr1_energy_smeared = [array.array('f', [0.0] * 256) for i in range(10)]
# tr2_n_hits = array.array('i', [0])
# tr2_energy_smeared = [array.array('f', [0.0] * 256) for i in range(10)]
cal_n_hits = array.array('i', [0])
cal_energy_smeared = [array.array('f', [0.0] * 256) for i in range(10)]


# t_output.Branch('tr1_n_hits', tr1_n_hits, 'tr1_n_hits/I')
# for _ in range(10):
#     t_output.Branch('tr1_energy_smeared{}'.format(_+1), tr1_energy_smeared[_], 'tr1_energy_smeared{}[tr1_n_hits]/F'.format(_+1))
#
# t_output.Branch('tr2_n_hits', tr2_n_hits, 'tr2_n_hits/I')
# for _ in range(10):
#     t_output.Branch('tr2_energy_smeared{}'.format(_+1), tr2_energy_smeared[_], 'tr2_energy_smeared{}[tr2_n_hits]/F'.format(_+1))

t_output.Branch('cal_n_hits', cal_n_hits, 'cal_n_hits/I')
for _ in range(6, 7):
    t_output.Branch('cal_energy_smeared{}'.format(_+1), cal_energy_smeared[_], 'cal_energy_smeared{}[cal_n_hits]/F'.format(_+1))


for idx, event in enumerate(t_input):
    if (idx % 1000) == 0:
        print("Event: ", idx)
    # if idx == 10000:
    #     break
    # tr1_n_hits[0] = event.tr1_n_hits
    # for i in range(tr1_n_hits[0]):
    #     sector = event.tr1_hit_sector[i]
    #     pad = event.tr1_hit_pad[i]
    #     layer = event.tr1_hit_layer[i]
    #     energy = event.tr1_hit_energy[i]
    #     for j in range(10):
    #         rand_noise = np.random.normal(0., 0.1*(j+1)*noise[sector][pad][layer])
    #         tr1_energy_smeared[j][i] = energy + rand_noise
    #
    # tr2_n_hits[0] = event.tr2_n_hits
    # for i in range(tr2_n_hits[0]):
    #     sector = event.tr2_hit_sector[i]
    #     pad = event.tr2_hit_pad[i]
    #     layer = event.tr2_hit_layer[i]
    #     energy = event.tr2_hit_energy[i]
    #     for j in range(10):
    #         rand_noise = np.random.normal(0., 0.1*(j+1)*noise[sector][pad][layer])
    #         tr2_energy_smeared[j][i] = energy + rand_noise

    cal_n_hits[0] = event.cal_n_hits
    for i in range(cal_n_hits[0]):
        sector = event.cal_hit_sector[i]
        pad = event.cal_hit_pad[i]
        layer = event.cal_hit_layer[i]
        energy = event.cal_hit_energy[i]
        for j in range(6,7):
            rand_noise = np.random.normal(0., 0.1*(j+1)*noise[sector][pad][layer])
            cal_energy_smeared[j][i] = energy + rand_noise

    t_output.Fill()

f_output.Write()
f_output.Close()
