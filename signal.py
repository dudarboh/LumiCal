import numpy as np
from apv_data import position, calib_energy


class Signal:
    '''
    Signal in the detector
    Attributes:
    sector, pad, layer - position of the signal
    energy - energy of the signal
    neighbor - the most energetic neighbor of the signal or itself.
    x, y - coordinates of signal in x, y.
    cluster - cluster index to which the signal was asigned

    position() - method to calculate position of the signal
    based on apv's id and channel where signal was recorded
    calib_energy - method to calculate calibrated energy of the
    signal.It uses: calibration files for this apv, recorded signal
    value of ADC (Voltage) in apv.
    '''

    def __init__(self, sector, pad, energy):
        # Calculate all parameters of the signal when object created
        self.sector = sector
        self.pad = pad
        self.energy = energy

        self.neighbor = -1
        self.cluster = -1
        # Transform pad/sector numbers into x,y coordinates
        rho = 80+0.9+1.8*self.pad
        phi = np.pi/2+np.pi/12-np.pi/48-np.pi/24*self.sector
        self.x = rho*np.cos(phi)
        self.y = rho*np.sin(phi)


def extract_signal(event):
    '''
    Write something clever here.
    '''

    # Read needed branches from the input ROOT file.
    id_arr = event.apv_id
    channel_arr = event.apv_ch
    signal_arr = event.apv_signal_maxfit
    apv_nn_output = event.apv_nn_output
    apv_fit_tau = event.apv_fit_tau
    apv_fit_t0 = event.apv_fit_t0
    apv_bint1 = event.apv_bint1

    # Lists for signals
    signals_calorimeter = []
    signals_tracker1 = []
    signals_tracker2 = []

    # Loop through all signals(hits) in the event
    for hit in range(len(id_arr)):
        if (apv_fit_tau[hit] < 1 or apv_fit_tau[hit] > 3
           or signal_arr[hit] > 2000.
           or apv_fit_t0[hit] < (apv_bint1[hit]-2.7)
           or apv_fit_t0[hit] > (apv_bint1[hit]-0.5)):
            continue

        # Calculate position
        sector, pad, layer = position(id_arr[hit], channel_arr[hit])

        if (sector == 0 or sector == 3 or layer == 7
           or (sector == 1 and pad < 20)
           or (sector == 2 and pad < 20)
           or sector < 0  # This one is changed due to python C++ difference in %.
           or (layer < 2 and (signal_arr[hit] < 0. or apv_nn_output[hit] < 0.5))):
            continue

        energy = calib_energy(id_arr[hit], signal_arr[hit])
        # Ignore noisy cells in calorimeter
        if layer >= 2 and (energy < 1.4 or apv_nn_output[hit] < 0.5):
            continue

        # Choose what is data_list(tracker1/2,calorimeter)
        if layer == 0:
            data_list = signals_tracker1
        elif layer == 1:
            data_list = signals_tracker2
        else:
            data_list = signals_calorimeter

        # If signal with this position already in the list: just add energy
        # Else: Add this signal to list and add energy
        for item in data_list:
            if (item.sector, item.pad) == (sector, pad):
                item.energy += energy
                break
        else:
            data_list.append(Signal(sector, pad, energy))

    # If no signals in calorimeter - skip event
    # if len(signals_calorimeter) == 0:
    #    return False

    # Sort signals by energy
    signals_calorimeter = sorted(signals_calorimeter, key=lambda x: x.energy, reverse=True)
    signals_tracker1 = sorted(signals_tracker1, key=lambda x: x.energy, reverse=True)
    signals_tracker2 = sorted(signals_tracker2, key=lambda x: x.energy, reverse=True)

    # If everything is ok
    return signals_tracker1, signals_tracker2, signals_calorimeter
