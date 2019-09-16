import numpy as np
from ROOT import gStyle, TGeoMaterial, TGeoMedium, TGeoTranslation, TGeoManager, TGeoRotation


def main():

    gStyle.SetCanvasPreferGL(True)
    geom = TGeoManager("LumiCal", "Lumical Sensor")

    matVacuum = TGeoMaterial("Vacuum", 0, 0, 0)
    matAl = TGeoMaterial("Al", 26.98, 13, 2.7)

    Vacuum = TGeoMedium("Vacuum", 1, matVacuum)
    Al = TGeoMedium("Root Material", 2, matAl)

    rotation = TGeoRotation('rot', 90., 0., 0)

    top = geom.MakeBox("TOP", Vacuum, 100., 100., 200.)
    geom.SetTopVolume(top)

    sensor = geom.MakeBox("Sensor", Vacuum, 90, 90, 190)
    sensor.SetVisibility(False)

    pad_vols = [[0. for _ in range(64)] for __ in range(4)]

    for sector in range(4):
        for pad in range(64):
            r_in_pad = 80. + 1.8 * pad
            r_out_pad = r_in_pad + 1.2
            phi_min = -15 + 7.5 * sector
            phi_max = phi_min + 7.

            pad_vols[sector][pad] = geom.MakeTubs("pad_vol_{}_{}".format(sector, pad), Al, r_in_pad, r_out_pad, 0.320, phi_min, phi_max)
            pad_vols[sector][pad].SetLineColor(2)

            sensor.AddNode(pad_vols[sector][pad], sector * 64 + pad, rotation)

    for _ in range(45, 52):
        pad_vols[1][_].SetLineColor(2)

    top.AddNode(sensor, 1, TGeoTranslation(0., 0., 0.))

    geom.CloseGeometry()
    top.Draw("ogl")
    input('wait')

main()
