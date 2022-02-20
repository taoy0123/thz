'''
Calculates signal power spectra in dB
Created by Yu Heng Tao 20 Feb 2022
'''
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import rfft, rfftfreq
from tkinter.filedialog import askopenfilenames, asksaveasfile
from pathlib import Path
import os

c = 299792458

class Signal:
    def __init__(self, amplitude, phase, frequency, **kwargs):
        self.amplitude = amplitude
        self.phase = phase
        self.frequency = frequency

def read_data(filenames):
    dfs = []
    for filename in filenames:
        _, file_extension = os.path.splitext(filename)
        if file_extension == '.txt':
            sep = '\t'
            skip = 7
        else:
            sep = ','
            skip = 3
        df = pd.read_csv(filename, sep=sep, skiprows=skip, usecols=[0,1], names=['ps', 'a.u.'], engine='python')
        dfs.append(df)
    return list(zip(dfs, filenames))

def create_signal(df):
    '''returns the fft of signal in amplitude and phase as a signal object'''
    
    t = df['ps'].values * 10**-12
    xt = df['a.u.'].values
    dT = t[2] - t[1]
    
    xf = rfft(xt)
    f = rfftfreq(len(xt), dT)
    
    amplitude = abs(xf)
    phase = np.unwrap(np.angle(xf))

    # shift all phase below zero
    ind = int(1*10**12 / f[1])
    gradient = (phase[ind+1]-phase[ind]) / (f[ind+1]-f[ind])
    y_int = phase[ind] - (gradient * f[ind])
    phase -= y_int

    return Signal(amplitude[1:], phase[1:], f[1:])

def calc_power(signal):
    '''returns the power spectrum in dB as an array'''
    return 20*np.log10(signal.amplitude/signal.amplitude.max())

def get_filenames():
    '''get names of csv files'''
    filenames = askopenfilenames()
    return filenames

if __name__ == '__main__':

    # read data
    filenames = get_filenames()
    measurements = read_data(filenames)

    # FFT
    thz_signals = []
    for measurement in measurements:
        thz_signal = create_signal(measurement[0])
        thz_signal.label =  Path(measurement[1]).stem
        thz_signals.append(thz_signal)
       
    # calculate power spectra
    powers=[]
    for thz_signal in thz_signals:
        thz_signal.power = calc_power(thz_signal)

    # plots for display
    fig, axs = plt.subplots()

    for thz_signal in thz_signals:
        axs.plot(thz_signal.frequency*10**-12, thz_signal.power, label=thz_signal.label)

    axs.set_xlim([0, 10])
    axs.set_xlabel('Frequency (THz)')
    axs.set_ylabel('Power (dB)')
    axs.legend()
    
    plt.tight_layout()
    plt.show()

    # save results
    export_list = []
    column_names = []
    for thz_signal in thz_signals:
        export_list.append(pd.DataFrame(thz_signal.frequency*10**-12))
        export_list.append(pd.DataFrame(thz_signal.power))
        column_names.append('THz')
        column_names.append(thz_signal.label)
    export_df = pd.concat(export_list, axis=1)
    export_df.columns = column_names

    try:
        with asksaveasfile(mode='w', defaultextension=".csv") as file:
            export_df.to_csv(file.name, index=False)

    except AttributeError:
        # if user cancels save, filedialog returns None rather than a file object, and the 'with' will raise an error
        print("The user cancelled save")