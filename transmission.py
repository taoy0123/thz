'''
Computes complex refractive index of sample given the sample and the reference time-domain electric fields
Created by Yu Heng Tao 29 Nov 2021
'''

import pandas as pd
import numpy as np
np.seterr(divide='ignore', invalid='ignore') #ignore divide by zero warning
import matplotlib.pyplot as plt
from scipy import signal
from scipy.fft import rfft, rfftfreq
import tkinter as tk
from tkinter import Tk, filedialog
from tkinter.filedialog import askopenfilename
from tkinter.messagebox import showinfo


c = 299792458

class Signal:
    def __init__(self, amplitude, phase, frequency, n=1, k=0, **kwargs):
        self.amplitude = amplitude
        self.phase = phase
        self.frequency = frequency
        self.n = n
        self.k = k

def read_data(filename, skip, sep):
    df = pd.read_csv(filename, sep=sep, skiprows=skip, names=['Delays (ps)', 'Signal (a.u.)'], engine='python')
    return df

def create_signal(df, window=1):
    '''returns the fft of signal in amplitude and phase as a signal object'''
    
    t = df['Delays (ps)'].values * 10**-12
    xt = df['Signal (a.u.)'].values
    
    T = max(t) - min(t)
    dT = T/(len(t)-1)
    
    
    xf = rfft(xt*window)
    f = rfftfreq(len(xt), dT)
    
    phase = np.angle(xf)
    
    return Signal(abs(xf), np.unwrap(phase), f)

def calc_power(signal):
    '''returns the power spectrum in dB as an array'''
    return 20*np.log10(signal.amplitude/signal.amplitude.max())

def calc_complex_n(sample, reference):
    '''returns the absorption and n of the sample as two arrays'''
    A = sample.amplitude / reference.amplitude
    phi = sample.phase - reference.phase

    w = 2*np.pi*sample.frequency

    d = sample.thickness
    n_ref = reference.n
    
    n_sample = 1 - c*phi/d/w
    
    t_co = 4*n_sample*n_ref/(n_sample+n_ref)**2
    k_sample = -c/d/w*np.log(A/t_co)
    
    return n_sample, k_sample

def get_filenames():
    '''get names of csv files'''
    Tk().withdraw()
    showinfo(message='Select sample file')  
    sample_filename = askopenfilename()

    showinfo(message='Select reference file')
    reference_filename = askopenfilename() 

    return sample_filename, reference_filename

def find_peak(df, start_ps=None):
    if start_ps == None:
        start_ps = df.iloc[0,0]
    peak_ind = df.loc[df['Delays (ps)'] >= start_ps,'Signal (a.u.)'].idxmax()
    peak_ps = df.iloc[peak_ind,0]
    return peak_ind, peak_ps

def shift_data(df, offset=0):
    shifted_df = df[:].copy()
    shifted_df.iloc[offset:,1] = df.iloc[:-offset,1]
    shifted_df.iloc[0:offset, 1] = 0
    return shifted_df

def create_window(df, second_peak_ps, offset=0):
    '''create time-window to exclude FP peaks'''
    # butterworth window 50th order, 3 ps by trial and error
    b, a = signal.butter(50, second_peak_ps-3-offset, analog=True)
    w, h = signal.freqs(b, a, worN=df['Delays (ps)'].values)
    if second_peak_ps < 0:
        butter_window = 1-abs(h)
    else:
        butter_window = abs(h)
    
    return butter_window
 
def get_user_input():
    global sample_thickness
    global reference_n
    global skip
    global sep
    sample_thickness = float(thickness_entry.get())*10**-6
    reference_n = float(reference_n_entry.get())
    skip = int(skip_entry.get())
    sep = sep_entry.get()
    root.quit()

if __name__ == '__main__':
    
    # get user inputs
    root = tk.Tk()
    var = tk.StringVar()
    var.set('teraview')

    sample_filename, reference_filename = get_filenames()

    tk.Label(root, text='Delimiter:  ').grid(row=0, column=0)
    sep_entry = tk.Entry(root)
    sep_entry.grid(row=0, column=1)

    tk.Label(root, text='No. header rows:  ').grid(row=1, column=0)
    skip_entry = tk.Entry(root)
    skip_entry.grid(row=1, column=1)

    tk.Label(root, text='Sample thickness (um):  ').grid(row=2, column=0)
    thickness_entry = tk.Entry(root)
    thickness_entry.grid(row=2, column=1)
    
    tk.Label(root, text='Reference refractive index:  ').grid(row=3, column=0)
    reference_n_entry = tk.Entry(root)
    reference_n_entry.grid(row=3, column=1)

    tk.Button(root, text='Calculate', command=get_user_input).grid(row=4, column=1)
    root.mainloop()
    root.withdraw()
    
    # read data
    df_sample = read_data(sample_filename, skip, sep)
    df_reference = read_data(reference_filename, skip, sep)
    
    # find peaks
    main_peak_ind_ref, main_peak_ps_ref = find_peak(df_reference)
    main_peak_ind_sample, main_peak_ps_sample = find_peak(df_sample)

    main_peak_delay_ind = main_peak_ind_sample - main_peak_ind_ref
    main_peak_delay_ps = main_peak_ps_sample - main_peak_ps_ref

    # shift reference to find FP
    df_shifted_ref = shift_data(df_reference, main_peak_delay_ind)
    # get covariance and find minimum
    df_diff = df_sample[:].copy()
    df_diff.iloc[:,1] = df_sample['Signal (a.u.)'].values - df_shifted_ref['Signal (a.u.)'].values
    
    # start searching for second peak 12 ps after main peak
    start_ps = main_peak_ps_sample + 5
    second_peak_ind_sample, second_peak_ps_sample = find_peak(df_diff, start_ps=start_ps)
    print(f'Second peak at {second_peak_ps_sample:.2f} ps')
    ax = df_sample.plot(x='Delays (ps)', y='Signal (a.u.)')
    df_shifted_ref.plot(x='Delays (ps)', y='Signal (a.u.)', ax=ax)
    plt.axvline(x=second_peak_ps_sample,color='r')
    plt.legend(['Sample', 'Shifted reference'])
    plt.show()

    # ask user for alternative FP peak location
    second_peak_ps_sample_corrected = input('Enter alternative FP peak location (ps):  ')
    if second_peak_ps_sample_corrected != '':
        second_peak_ps_sample = float(second_peak_ps_sample_corrected)
    
    sample_window = create_window(df_sample ,second_peak_ps_sample)
    ref_window = create_window(df_sample, second_peak_ps_sample, offset=main_peak_delay_ps)
    
    sample = create_signal(df_sample, window=sample_window)
    reference = create_signal(df_reference, window=ref_window)
    
    sample.thickness = sample_thickness
    reference.n = reference_n
    
    # calculate optical parameters
    sample.n, sample.k = calc_complex_n(sample, reference)
    sample.absorption = 2*(2*np.pi*sample.frequency)/c*sample.k/100
    
    # TODO: calculate n obtained from first FP
    # TODO: calculate estimated sample thickness from the difference between n
    
    # calculate power spectra
    sample.power = calc_power(sample)
    reference.power = calc_power(reference)
    
    # plots for display
    fig, axs = plt.subplots(2, 2)
    
    axs[0,0].plot(sample.frequency*10**-12, sample.n)
    axs[1,0].plot(sample.frequency*10**-12, sample.absorption)

    axs[0,1].plot(df_sample['Delays (ps)'].values, df_sample['Signal (a.u.)'].values, label='Sample')
    axs[0,1].plot(df_shifted_ref['Delays (ps)'].values, df_shifted_ref['Signal (a.u.)'].values, label='Shifted reference')
    axs[0,1].plot(df_sample['Delays (ps)'].values, sample_window, label='Sample window')
    axs[0,1].axvline(x=second_peak_ps_sample, color='r')

    axs[1,1].plot(sample.frequency*10**-12, sample.power, label='Sample')
    axs[1,1].plot(reference.frequency*10**-12, reference.power, label='Reference') 
    
    axs[0,0].set_xlim([0.2, 5])
    axs[0,0].set_ylim([0,3])
    axs[1,0].set_xlim([0.2, 5])
    axs[1,0].set_ylim(bottom=0)
    axs[1,1].set_xlim([0, 10])
    
    axs[0,0].set_xlabel('Frequency (THz)')
    axs[0,0].set_ylabel('Refractive index')
    axs[1,0].set_xlabel('Frequency (THz)')
    axs[1,0].set_ylabel('Absorption coefficient ($cm^{-1}$)')
    axs[0,1].set_xlabel('Delays (ps)')
    axs[0,1].set_ylabel('Signal (a.u.)')
    axs[1,1].set_xlabel('Frequency (THz)')
    axs[1,1].set_ylabel('Power (dB)')
    
    axs[0,1].legend()
    axs[1,1].legend()
    
    plt.tight_layout()
    plt.show()

    # save results
    export_array = np.stack((sample.frequency, sample.n, sample.absorption)).T[1:]
    export_df = pd.DataFrame(export_array, columns=['f','n','a'])

    try:
        with filedialog.asksaveasfile(mode='w', defaultextension=".csv") as file:
            export_df.to_csv(file.name, index=False)

    except AttributeError:
        # if user cancels save, filedialog returns None rather than a file object, and the 'with' will raise an error
        print("The user cancelled save")