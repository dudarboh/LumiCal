from ROOT import TMath
from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.morphology import generate_binary_structure, binary_erosion


# ###

# ## coverts APV's channel to it's pad position (id) by the following scheme:
# Sectors:
#      0           1         2         3
#  _________________________________________
# |   63     |    127   |   191   |   255   |
# |-----------------------------------------|
#  |   62    |    126   |   190   |  254   |
#  |---------------------------------------|
#   |   61    |   125   |   189  |  253   |
#   |-------------------------------------|
#
#          .........................
#          .........................
#
#          |  0  |  64 | 128 | 192 |
#          |_____|_____|_____|_____|

# Number of APVs used in experiment
# n_apvs = 16
# Number of channels of one APV
n_channels = 128
# Number of layers in experiment
n_layers = 6
# Number of sectors of the layer in the experiment
n_sectors = 4
# Number of pads of one sector
n_pads = 64

# Inner radius of detector and length of 1 pad in radial direction
r_inner = 80  # mm
r_pad = 1.8  # mm

apv_maps = {}


# tb15_master
apv_maps['tb15_master'] = []
for i in range(n_channels):
    if i < n_channels-1:
        apv_maps['tb15_master'] += [190-i if i < 63 else i+129]
    else:
        apv_maps['tb15_master'] += [-1]

# tb15_slave
apv_maps['tb15_slave'] = [-1, 62, 63, 60, 61, 58, 59, 56, 57, 54, 55, 52, 53,
                          50, 51, 48, 49, 46, 47, 44, 45, 42, 43, 40, 41, 38,
                          39, 36, 37, 34, 35, 32, 33, 30, 31, 28, 29, 26, 27,
                          24, 25, 22, 23, 20, 21, 18, 19, 16, 17, 14, 15, 12,
                          13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1, 65, 64, 67,
                          66, 69, 68, 71, 70, 73, 72, 75, 74, 77, 76, 79, 78,
                          81, 80, 83, 82, 85, 84, 87, 86, 89, 88, 91, 90, 93,
                          92, 95, 94, 97, 96, 99, 98, 101, 100, 103, 102, 105,
                          104, 107, 106, 109, 108, 111, 110, 113, 112, 115,
                          114, 117, 116, 119, 118, 121, 120, 123, 122, 125,
                          124, 127]

# tb16_master_tab_divider
apv_maps['tb16_master_divider'] = [-1, 255, 254, 253, 252, 251, 250, 249, 248,
                                   247, 246, 245, 244, 243, 242, 241, 240, 239,
                                   238, 237, 236, 235, 234, 233, 232, 231, 230,
                                   229, 228, 227, 226, 225, 224, 223, 222, 221,
                                   220, 219, 218, 217, 216, 215, 214, 213, 212,
                                   211, 210, 209, 208, 207, 206, 205, 204, 203,
                                   202, 201, 200, 199, 198, 197, 196, 195, 194,
                                   193, 192, 128, 129, 130, 131, 132, 133, 134,
                                   135, 136, 137, 138, 139, 140, 141, 142, 143,
                                   144, 145, 146, 147, 148, 149, 150, 151, 152,
                                   153, 154, 155, 156, 157, 158, 159, 160, 161,
                                   162, 163, 164, 165, 166, 167, 168, 169, 170,
                                   171, 172, 173, 174, 175, 176, 177, 178, 179,
                                   180, 181, 182, 183, 184, 185, 186, 187, 188,
                                   189, 191]

# tb16_slave_divider
apv_maps['tb16_slave_divider'] = [126, 124, 125, 122, 123, 120, 121, 118, 119,
                                  116, 117, 114, 115, 112, 113, 110, 111, 108,
                                  109, 106, 107, 104, 105, 102, 103, 100, 101,
                                  98, 99, 96, 97, 94, 95, 92, 93, 90, 91, 88,
                                  89, 86, 87, 84, 85, 82, 83, 80, 81, 78, 79,
                                  76, 77, 74, 75, 72, 73, 70, 71, 68, 69, 66,
                                  67, 64, 65, 1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11,
                                  10, 13, 12, 15, 14, 17, 16, 19, 18, 21, 20,
                                  23, 22, 25, 24, 27, 26, 29, 28, 31, 30, 33,
                                  32, 35, 34, 37, 36, 39, 38, 41, 40, 43, 42,
                                  45, 44, 47, 46, 49, 48, 51, 50, 53, 52, 55,
                                  54, 57, 56, 59, 58, 61, 60, 63, 62, -1]

# tb16_master_tab_divider
apv_maps['tb16_master_tab_divider'] = [191, 189, 188, 187, 186, 185, 184, 183,
                                       182, 181, 180, 179, 178, 177, 176, 175,
                                       174, 173, 172, 171, 170, 169, 168, 167,
                                       166, 165, 164, 163, 162, 161, 160, 159,
                                       158, 157, 156, 155, 154, 153, 152, 151,
                                       150, 149, 148, 147, 146, 145, 144, 143,
                                       142, 141, 140, 139, 138, 137, 136, 135,
                                       134, 133, 132, 131, 130, 129, 128, 192,
                                       193, 194, 195, 196, 197, 198, 199, 200,
                                       201, 202, 203, 204, 205, 206, 207, 208,
                                       209, 210, 211, 212, 213, 214, 215, 216,
                                       217, 218, 219, 220, 221, 222, 223, 224,
                                       225, 226, 227, 228, 229, 230, 231, 232,
                                       233, 234, 235, 236, 237, 238, 239, 240,
                                       241, 242, 243, 244, 245, 246, 247, 248,
                                       249, 250, 251, 252, 253, 254, 255, -1]

# tb16_slave_tab_divider
apv_maps['tb16_slave_tab_divider'] = [126, 124, 125, 122, 123, 120, 121, 118,
                                      119, 116, 117, 114, 115, 112, 113, 110,
                                      111, 108, 109, 106, 107, 104, 105, 102,
                                      103, 100, 101, 98, 99, 96, 97, 94, 95,
                                      92, 93, 90, 91, 88, 89, 86, 87, 84, 85,
                                      82, 83, 80, 81, 78, 79, 76, 77, 74, 75,
                                      72, 73, 70, 71, 68, 69, 66, 67, 64, 65,
                                      1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13,
                                      12, 15, 14, 17, 16, 19, 18, 21, 20, 23,
                                      22, 25, 24, 27, 26, 29, 28, 31, 30, 33,
                                      32, 35, 34, 37, 36, 39, 38, 41, 40, 43,
                                      42, 45, 44, 47, 46, 49, 48, 51, 50, 53,
                                      52, 55, 54, 57, 56, 59, 58, 61, 60, 63,
                                      62, -1]

# ### Calibration files for energy in every APV
# path = '/data/alzta/aborysov/tb_2016_data/code/lumical_clust/fcalib/'
calib_path = '../calibration/'
calib_file_names = ['calibration_apv_0.txt', 'calibration_apv_1.txt',
                    'calibration_apv_2.txt',  'calibration_apv_3.txt',
                    'calibration_apv_4.txt',  'calibration_apv_5.txt',
                    'calibration_apv_6.txt',  'calibration_apv_7.txt',
                    'calibration_apv_8.txt',  'calibration_apv_9.txt',
                    'calibration_apv_10.txt', 'calibration_apv_11.txt',
                    'calibration_apv_12.txt', 'calibration_apv_13.txt',
                    'calibration_apv_14.txt', 'calibration_apv_15.txt']


def langaufun(x, par):
    # Fit parameters:
    # par[0]=Width (scale) parameter of Landau density
    # par[1]=Most Probable (MP, location) parameter of Landau density
    # par[2]=Total area (integral -inf to inf, normalization constant)
    # par[3]=Width (sigma) of convoluted Gaussian function
    #
    # In the Landau distribution (represented by the CERNLIB approximation),
    # the maximum is located at x=-0.22278298 with the location parameter=0.
    # This shift is corrected within this function, so that the actual
    # maximum is identical to the MP parameter.

    # Numeric constants
    invsq2pi = 0.3989422804014  # (2 pi)^(-1/2)
    mpshift = -0.22278298  # Landau maximum location

    # Control constants
    np = 1000  # number of convolution steps
    sc = 5.  # convolution extends to +-sc Gaussian sigmas

    # Variables
    summ = 0

    # MP shift correction
    mpc = par[1]-mpshift*par[0]

    # Range of convolution integral
    xlow = x[0]-sc*par[3]
    xupp = x[0]+sc*par[3]

    step = (xupp-xlow)/np

    # Convolution integral of Landau and Gaussian by sum
    # Fix divide by zero error
    if par[0] == 0:
        par[0] = 1e-6
    if par[3] == 0:
        par[3] = 1e-6
    for i in range(np//2):
        xx = xlow+(i+0.5)*step
        fland = TMath.Landau(xx, mpc, par[0])/par[0]
        summ += fland*TMath.Gaus(x[0], xx, par[3])

        xx = xupp-(i+0.5)*step
        fland = TMath.Landau(xx, mpc, par[0])/par[0]
        summ += fland*TMath.Gaus(x[0], xx, par[3])

    return par[2]*step*summ*invsq2pi/par[3]


def detect_peaks(array):
    """
    Takes an 2d array and detect the peaks using the local maximum filter.
    Returns a boolean mask of the peaks (i.e. 1 when
    the pixel's value is the neighborhood maximum, 0 otherwise)
    """

    # define an 8-connected neighborhood - count diagonal elements as
    # neighbors
    neighborhood = generate_binary_structure(2, 2)

    # apply the local maximum filter; all pixel of maximal value
    # in their neighborhood are set to 1
    local_max = maximum_filter(array, footprint=neighborhood) == array
    # local_max is a mask that contains the peaks we are
    # looking for, but also the background.
    # In order to isolate the peaks we must remove
    # the background from the mask.

    # we create the mask of the background
    background = (array == 0)

    # a little technicality: we must erode the background in order to
    # successfully subtract it form local_max, otherwise a line will
    # appear along the background border
    # (artifact of the local maximum filter)
    eroded_background = binary_erosion(background, structure=neighborhood,
                                       border_value=1)

    # we obtain the final mask, containing only peaks,
    # by removing the background from the local_max mask (xor operation)
    result = local_max ^ eroded_background

    return result
