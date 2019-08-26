from ROOT import TFile, TGraphErrors, TH2F, gStyle
import numpy as np
from itertools import islice


class ApvMaps:
    ''' Maps of apv channels to pad number. Are taken from Sasha  analysis code'''
    tb15_master = [190 - i if i < 63 else i + 129 for i in range(127)] + [-1]

    tb15_slave = [-1, 62, 63, 60, 61, 58, 59, 56, 57, 54, 55, 52, 53, 50, 51,
                  48, 49, 46, 47, 44, 45, 42, 43, 40, 41, 38, 39, 36, 37, 34,
                  35, 32, 33, 30, 31, 28, 29, 26, 27, 24, 25, 22, 23, 20, 21,
                  18, 19, 16, 17, 14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5,
                  2, 3, 0, 1, 65, 64, 67, 66, 69, 68, 71, 70, 73, 72, 75, 74,
                  77, 76, 79, 78, 81, 80, 83, 82, 85, 84, 87, 86, 89, 88, 91,
                  90, 93, 92, 95, 94, 97, 96, 99, 98, 101, 100, 103, 102, 105,
                  104, 107, 106, 109, 108, 111, 110, 113, 112, 115, 114, 117,
                  116, 119, 118, 121, 120, 123, 122, 125, 124, 127]

    tb15_slave = [-1, 62, 63, 60, 61, 58, 59, 56, 57, 54, 55, 52, 53, 50, 51,
                  48, 49, 46, 47, 44, 45, 42, 43, 40, 41, 38, 39, 36, 37, 34,
                  35, 32, 33, 30, 31, 28, 29, 26, 27, 24, 25, 22, 23, 20, 21,
                  18, 19, 16, 17, 14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2,
                  3, 0, 1, 65, 64, 67, 66, 69, 68, 71, 70, 73, 72, 75, 74, 77,
                  76, 79, 78, 81, 80, 83, 82, 85, 84, 87, 86, 89, 88, 91, 90, 93,
                  92, 95, 94, 97, 96, 99, 98, 101, 100, 103, 102, 105, 104, 107,
                  106, 109, 108, 111, 110, 113, 112, 115, 114, 117, 116, 119, 118,
                  121, 120, 123, 122, 125, 124, 127]

    tb16_master_divider = [-1, 255, 254, 253, 252, 251, 250, 249, 248, 247, 246,
                           245, 244, 243, 242, 241, 240, 239, 238, 237, 236, 235,
                           234, 233, 232, 231, 230, 229, 228, 227, 226, 225, 224,
                           223, 222, 221, 220, 219, 218, 217, 216, 215, 214, 213,
                           212, 211, 210, 209, 208, 207, 206, 205, 204, 203, 202,
                           201, 200, 199, 198, 197, 196, 195, 194, 193, 192, 128,
                           129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139,
                           140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150,
                           151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161,
                           162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172,
                           173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183,
                           184, 185, 186, 187, 188, 189, 191]

    tb16_slave_divider = [126, 124, 125, 122, 123, 120, 121, 118, 119, 116, 117,
                          114, 115, 112, 113, 110, 111, 108, 109, 106, 107, 104,
                          105, 102, 103, 100, 101, 98, 99, 96, 97, 94, 95, 92, 93,
                          90, 91, 88, 89, 86, 87, 84, 85, 82, 83, 80, 81, 78, 79, 76,
                          77, 74, 75, 72, 73, 70, 71, 68, 69, 66, 67, 64, 65, 1, 0, 3,
                          2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14, 17, 16, 19, 18,
                          21, 20, 23, 22, 25, 24, 27, 26, 29, 28, 31, 30, 33, 32, 35,
                          34, 37, 36, 39, 38, 41, 40, 43, 42, 45, 44, 47, 46, 49, 48,
                          51, 50, 53, 52, 55, 54, 57, 56, 59, 58, 61, 60, 63, 62, -1]

    tb16_master_tab_divider = [191, 189, 188, 187, 186, 185, 184, 183, 182, 181, 180,
                               179, 178, 177, 176, 175, 174, 173, 172, 171, 170, 169,
                               168, 167, 166, 165, 164, 163, 162, 161, 160, 159, 158,
                               157, 156, 155, 154, 153, 152, 151, 150, 149, 148, 147,
                               146, 145, 144, 143, 142, 141, 140, 139, 138, 137, 136,
                               135, 134, 133, 132, 131, 130, 129, 128, 192, 193, 194,
                               195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205,
                               206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216,
                               217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227,
                               228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238,
                               239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249,
                               250, 251, 252, 253, 254, 255, -1]

    tb16_slave_tab_divider = [126, 124, 125, 122, 123, 120, 121, 118, 119, 116, 117,
                              114, 115, 112, 113, 110, 111, 108, 109, 106, 107, 104,
                              105, 102, 103, 100, 101, 98, 99, 96, 97, 94, 95, 92, 93,
                              90, 91, 88, 89, 86, 87, 84, 85, 82, 83, 80, 81, 78, 79,
                              76, 77, 74, 75, 72, 73, 70, 71, 68, 69, 66, 67, 64, 65,
                              1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14,
                              17, 16, 19, 18, 21, 20, 23, 22, 25, 24, 27, 26, 29, 28,
                              31, 30, 33, 32, 35, 34, 37, 36, 39, 38, 41, 40, 43, 42,
                              45, 44, 47, 46, 49, 48, 51, 50, 53, 52, 55, 54, 57, 56,
                              59, 58, 61, 60, 63, 62, -1]


class CalibGraphs:
    '''APV calibration. To convert Volts to MIPs'''
    # path on alzt.tau.ac.il server = '/data/alzta/aborysov/tb_2016_data/code/lumical_clust/fcalib/'
    calib_graphs = []

    @classmethod
    def get_calib_graphs(cls):
        path = "../calibration_files/"
        for i in range(16):
            calib_file = path + "calibration_apv_{}".format(i) + ".txt"

            # 1st point
            x = [0.]
            y = [0.]
            x_err = [1.e-5]
            y_err = [1.e-5]

            # Calibration x-y data is inverted
            with open(calib_file, 'r') as file:
                for line in islice(file, 1, None):
                    x.append(float(line.split('  ')[1]))
                    y.append(float(line.split('  ')[0]))
                    x_err.append(float(line.split('  ')[3]))
                    y_err.append(float(line.split('  ')[2]))

            x = np.array(x)
            y = np.array(y) * 19.206
            x_err = np.array(x_err)
            y_err = np.array(y_err) * 19.206

            cls.calib_graphs.append(TGraphErrors(len(x), x, y, x_err, y_err))


def position(apv_id, apv_channel):
    if apv_id < 4:
        apv_map = ApvMaps.tb15_slave if apv_id % 2 == 1 else ApvMaps.tb15_master

    elif apv_id >= 4 and apv_id < 14:
        apv_map = ApvMaps.tb16_slave_divider if apv_id % 2 == 1 else ApvMaps.tb16_master_divider

    elif apv_id == 14:
        apv_map = ApvMaps.tb16_master_tab_divider

    elif apv_id == 15:
        apv_map = ApvMaps.tb16_slave_tab_divider

    layer = apv_id // 2
    sector = apv_map[apv_channel] // 64
    pad = apv_map[apv_channel] % 64
    return sector, pad, layer


def calib_energy(apv_id, apv_signal):
    signal = apv_signal if apv_signal < 1450. else 1450.
    return CalibGraphs.calib_graphs[apv_id].Eval(signal)


file = TFile.Open("/home/FoxWise/Downloads/run741.root")
tree = file.pedestals

output_filie = TFile.Open("noise_histos.root", "RECREATE")
output_filie.cd()

CalibGraphs.get_calib_graphs()

sector_list = []
pad_list = []
layer_list = []
noise_list = []

histos = [TH2F("layer_{}".format(i), "Noise_in_LAYER_{}".format(i), 2, 1, 3, 44, 20, 64) for i in range(8)]
for h in histos:
    h.GetYaxis().SetTitle("pad")
    h.GetXaxis().SetTitle("sector")
for event in tree:
    for i in range(2048):
        sector, pad, layer = position(event.apv_id[i], event.apv_ch[i])
        noise = calib_energy(event.apv_id[i], event.apv_pedstd[i]) * 0.0885  # This is in MeV
        sector_list.append(sector)
        pad_list.append(pad)
        layer_list.append(layer)
        noise_list.append(noise)
        if (sector == 1 or sector == 2) and pad>20:
            histos[layer].Fill(sector, pad, event.apv_pedstd[i])

gStyle.SetOptStat(0)
for idx, h in enumerate(histos):
    h.Write("noise_layer_{}".format(idx))
input("wait")
with open('noise.txt', 'w+') as f:
    f.write("sector  pad  layer noise_in_MeV\n")
    for i in range(2048):
        f.write("%s  %s  %s %s\n" % (sector_list[i], pad_list[i], layer_list[i], noise_list[i]))
