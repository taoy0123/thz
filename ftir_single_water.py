import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

from tkinter.filedialog import askopenfilename, askopenfilenames
from tkinter.messagebox import showinfo
from tkinter import Tk, filedialog

from scipy.signal import find_peaks, butter, freqs
from scipy.fft import rfft, rfftfreq, irfft

def read_data(filename, sep='\t', header=None):
    '''reads data into dataframe'''
    if header is None:
        df = pd.read_csv(filename, sep=sep, names=['wavenumber', 'intensity'], engine='python')
    else:
        df = pd.read_csv(filename, sep=sep, header=header, engine='python')
    df = df.sort_values(df.columns[0])
    return df

def get_filename():
    '''get names of input files'''
    showinfo(message='Select sample file')  
    sample_filename = askopenfilename()
    # showinfo(message='Select reference file')
    # select lowest temperature water measured as reference
    reference_filename = os.path.join(os.path.dirname(__file__), '13_water_THz_44C_MM8_Globar_2cm-1_128Scans_ev-ev.0.txt')
    return sample_filename, reference_filename

def get_filenames():
    '''get names of csv files'''
    showinfo(message='Select sample files') 
    sample_filenames = askopenfilenames()
    showinfo(message='Select reference files')
    reference_filenames = askopenfilenames()
    return sorted(sample_filenames), sorted(reference_filenames)

def create_signal_list(filenames):
    dfs = []
    temperatures = []
    for filename in filenames:
        df = read_data(filename)
        dfs.append(df)
        index = filename.find('C')
        temperature = filename[index-2:index+1]
        temperatures.append(temperature)
    return list(zip(dfs, temperatures))

def filter_etalon(f, xf):
    '''filter out etalon oscillations in xf'''
    # convert nan to zero. nan due to -ve intensity
    np.nan_to_num(xf,copy=False)

    # exclude >400cm-1 frequencies to see etalon more clearly
    # trim_xf = xf[:385]
    # trim_f = f[:385]
    trim_xf = xf[:644]
    trim_f = f[:644]

    # fourier transform
    dF = f[2] - f[1]
    xt = rfft(trim_xf)
    t = rfftfreq(len(trim_xf), dF)

    # find peaks
    peaks, _ = find_peaks(abs(xt), prominence=0.2)

    # determine cutoff freq for bandstop filter
    tolerance = 0.01
    focused_peaks = []
    for i in peaks:
        if t[i] > 0.2 and t[i] < 0.3:
            focused_peaks.append(t[i])
    # if there is peak in the range
    if focused_peaks:
        cutoff_freqs = [focused_peaks[0]-tolerance, focused_peaks[-1]+tolerance]
        # create bandstop filter
        b, a = butter(7, cutoff_freqs, analog=True, btype='bandstop')
        # w, h = freqs(b, a, worN=rfftfreq(len(xf), dF))
        w, h = freqs(b, a, worN=t)
        window = abs(h)
    else:
        print('No peaks found.')
        window = np.ones(len(t))

    # apply filter to input signal
    filtered_xt = window*xt

    # plot results
    # plt.plot(t, abs(xt))
    # plt.plot(t[peaks], abs(xt)[peaks], 'x')
    # plt.plot(t, window)
    # plt.plot(t, abs(filtered_xt))
    # plt.legend(['Unfiltered', 'Detected Peaks', 'Filter Window', 'Filtered'])
    # plt.show()

    # inverse fourier transform the filtered signal
    filtered_xf = irfft(filtered_xt, len(trim_xf))

    return trim_f, filtered_xf

def vapor_correction(f, xf, vapor_spectrum):
    # locate the base of 150cm-1 vapor peak
    vapor_zero_index = abs(f - 147.55).argmin()

    # locate the top of 150cm-1 vapor peak
    vapor_max_index = abs(f - 150.44).argmin()

    # calculate the 150 cm-1 peak height in the signal
    vapor_amount = abs(xf[vapor_max_index] - xf[vapor_zero_index])

    # determine the amount of vapor in the signal
    scaling_factor = vapor_amount / 30.42173
    residual_vapor = scaling_factor * vapor_spectrum[:len(xf)]

    # subtract the residual vapor spectrum from the signal
    if xf[vapor_max_index] < xf[vapor_zero_index]:
        corrected_spectrum = xf + residual_vapor
    else:
        corrected_spectrum = xf - residual_vapor


    # plot results
    # plt.plot(f, xf)
    # plt.plot(f, residual_vapor)
    # plt.plot(f, corrected_spectrum)
    # plt.legend(['With vapor', 'Residual vapor', 'Without vapor'])
    # plt.show()

    return residual_vapor, corrected_spectrum

def calc_molar_a_water(delta_A, ref_a, density, M=18.0153, sample_thickness=0.0013):
    # calculate sample absorption coefficient
    delta_a = delta_A / sample_thickness
    sample_a = delta_a + ref_a

    molar_a_sample = sample_a / density /1000 * M

    return delta_a, sample_a, molar_a_sample

if __name__ == '__main__':
    # read in data
    sample_filename, reference_filename = get_filename()
    sample_df = read_data(sample_filename)
    reference_df = read_data(reference_filename)

    # read in water vapor data
    # data above 400 cm-1 fluctuates wildly and may be excluded
    target_path = os.path.join(os.path.dirname(__file__), 'WaterVapor_400.csv')
    vapor_df = read_data(target_path, sep=',')    

    frequency = sample_df['wavenumber'].values
    sample_intensity = sample_df['intensity'].values
    reference_intensity = reference_df['intensity'].values

    vapor_frequency = vapor_df['wavenumber'].values
    vapor_intensity = vapor_df['intensity'].values

    # compute difference absorbance
    delta_A = -np.log(sample_intensity/reference_intensity)

    # filter etalon
    # filtered_frequency, delta_A_filtered = filter_etalon(frequency, delta_A)
    # Bypass etalon filtering 
    filtered_frequency = frequency
    delta_A_filtered = delta_A

    # filter water vapor
    residual_vapor, delta_A_corrected = vapor_correction(filtered_frequency, delta_A_filtered, vapor_intensity)

    # plot difference absorbance
    # plt.plot(frequency, delta_A)
    # plt.plot(filtered_frequency, delta_A_corrected)
    # ax = plt.gca()
    # ax.set_ylim([-0.5,0.5])
    # ax.set_xlabel('Wavenumbers (cm$^{-1}$)')
    # ax.set_ylabel(r'$\Delta$ A')
    # plt.legend(['Raw', 'Corrected'])
    # plt.show()

    # read in water absorption and extinction data
    target_path = os.path.join(os.path.dirname(__file__), 'waterAbsorption_44.csv')
    water_a_df = read_data(target_path, sep=',', header=0)
    ref_a = water_a_df['44C']


    # Load density of the water sample based on its temperature
    # read in water densities
    target_path = os.path.join(os.path.dirname(__file__), 'WaterDensities.csv')
    water_densities = read_data(target_path, sep=',', header=0)
    # select density of water sample
    index = sample_filename[1:].rfind('C')
    sample_temperature = int(sample_filename[index-1:index+1].lstrip('_'))
    density = water_densities.loc[water_densities['degC'] == sample_temperature, 'g/cm3'].item()

    # # ask user input for the density of the water sample
    # density = float(input('Input sample density (g/cm3):  '))

    
    # calculate extinction coefficient
    delta_a_water, absorption_water, epsilon_water = calc_molar_a_water(delta_A_corrected, ref_a, density)

    # save results
    export_array = np.stack((filtered_frequency, delta_a_water, absorption_water, epsilon_water)).T[1:]
    export_df = pd.DataFrame(export_array, columns=['cm-1', 'delta a', 'water absorption', 'water epsilon'])

    try:
        with filedialog.asksaveasfile(mode='w', defaultextension=".csv") as file:
            export_df.to_csv(file.name, index=False)

    except TypeError:
        # if user cancels save, filedialog returns None rather than a file object, and the 'with' will raise an error
        print("The user cancelled save")