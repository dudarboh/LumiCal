import numpy as np
from ROOT import gStyle, TGeoMaterial, TGeoMedium, TGeoTranslation, TGeoManager, TGeoRotation



def bad_pad(sector, pad, layer):
    """Return true if pad is bad"""
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


def main(layer):

    # gStyle.SetCanvasPreferGL(True)
    geom = TGeoManager("LumiCal", "Lumical Sensor")

    matVacuum = TGeoMaterial("Vacuum", 0, 0, 0)
    matSi = TGeoMaterial("Si", 26.98, 13, 2.7)

    Vacuum = TGeoMedium("Vacuum", 1, matVacuum)
    Si = TGeoMedium("Root Material", 2, matSi)

    rotation = TGeoRotation('rot', 90., 0., 0)

    top = geom.MakeBox("TOP", Vacuum, 100., 100., 200.)
    geom.SetTopVolume(top)

    sensor = geom.MakeBox("Sensor", Si, 90, 90, 190)
    sensor.SetVisibility(False)

    pad_vols = [[0. for _ in range(64)] for __ in range(4)]

    for sector in range(4):
        for pad in range(64):
            r_in_pad = 80. + 1.8 * pad
            r_out_pad = r_in_pad + 1.2
            phi_min = -15 + 7.5 * sector
            phi_max = phi_min + 7.

            pad_vols[sector][pad] = geom.MakeTubs("pad_vol_{}_{}".format(sector, pad), Si, r_in_pad, r_out_pad, 0.320, phi_min, phi_max)
            if layer == 7:
                pad_vols[sector][pad].SetLineColor(2)
            elif pad < 20 or sector == 0 or sector == 3 or bad_pad(sector, pad, layer):
                pad_vols[sector][pad].SetLineColor(2)
            else:
                pad_vols[sector][pad].SetLineColor(11)
            sensor.AddNode(pad_vols[sector][pad], sector * 64 + pad, rotation)

    # for _ in range(45, 52):
    #     pad_vols[1][_].SetLineColor(3)

    top.AddNode(sensor, 1, TGeoTranslation(0., 0., 0.))

    geom.CloseGeometry()
    top.Draw("ogl")
    input('wait')


def main_full():

    # gStyle.SetCanvasPreferGL(True)
    geom = TGeoManager("LumiCal", "Lumical Sensor")

    matVacuum = TGeoMaterial("Vacuum", 0, 0, 0)
    matSi = TGeoMaterial("Si", 26.98, 13, 2.7)

    Vacuum = TGeoMedium("Vacuum", 1, matVacuum)
    Si = TGeoMedium("Root Material", 2, matSi)

    rotation = TGeoRotation('rot', 90., 0., 0)

    top = geom.MakeBox("TOP", Vacuum, 100., 100., 1000.)
    geom.SetTopVolume(top)

    sensor = geom.MakeBox("Sensor", Si, 90, 90, 190)
    sensor.SetVisibility(False)

    pad_vols = [[0. for _ in range(64)] for __ in range(4)]

    for sector in range(4):
        for pad in range(64):
            r_in_pad = 80. + 1.8 * pad
            r_out_pad = r_in_pad + 1.2
            phi_min = -15 + 7.5 * sector
            phi_max = phi_min + 7.

            pad_vols[sector][pad] = geom.MakeTubs("pad_vol_{}_{}".format(sector, pad), Si, r_in_pad, r_out_pad, 4., phi_min, phi_max)
            pad_vols[sector][pad].SetLineColor(11)
            if pad == 45 and sector == 3:
                pad_vols[sector][pad].SetLineColor(3)
            sensor.AddNode(pad_vols[sector][pad], sector * 64 + pad, rotation)

    # for _ in range(45, 52):
    #     pad_vols[1][_].SetLineColor(3)

    for _ in range(6):
        top.AddNode(sensor, _, TGeoTranslation(0., 0., _ * 10.))

    geom.CloseGeometry()
    top.Draw("ogl")
    input('wait')



main_full()
