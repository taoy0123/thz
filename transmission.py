'''
Computes complex refractive index of sample given the sample and the reference time-domain electric fields
Created by Yu Heng Tao 29 Nov 2021
'''

from cProfile import label
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
    df = pd.read_csv(filename, sep=sep, skiprows=skip, usecols=[0,1], names=['Delays (ps)', 'Signal (a.u.)'], engine='python')
    return df


def create_signal(df, window=1):
    '''returns the fft of signal in amplitude and phase as a signal object'''
    
    t = df['Delays (ps)'].values * 10**-12
    xt = df['Signal (a.u.)'].values
    dT = t[2] - t[1]
    
    xf = rfft(xt*window)
    f = rfftfreq(len(xt), dT)
    
    amplitude = abs(xf)
    angle = np.angle(xf)

    phase = np.unwrap(angle)

    # shift all phase below zero
    ind = int(1*10**12 / f[1])

    gradient = (phase[ind+1]-phase[ind]) / (f[ind+1]-f[ind])
    y_int = phase[ind] - (gradient * f[ind])
    phase -= y_int

    return Signal(amplitude[1:], phase[1:], f[1:])

def calc_power(signal, max):
    '''returns the power spectrum in dB as an array'''
    return 20*np.log10(signal.amplitude/max)


def find_ind(array, threshold):
    '''returns the index at which the threshold is first reached'''
    return np.where(array>=threshold)[0][0]

def calc_DR(signal):
    ''' returns the freq-domain dynamic range as an array'''
    # locate the start of noise floor
    # noise_ind_start = np.where(signal.frequency>=6*10**12)[0][0]
    noise_ind_start = find_ind(signal.frequency, 6*10**12)
    # locate the end of noise floor as some systems may have artefacts at very high frequencies
    try:
        # noise_ind_end = np.where(signal.frequency>15*10**12)[0][0] 
        noise_ind_end = find_ind(signal.frequency, 15*10**12)
    except IndexError:
        noise_floor = np.mean(signal.amplitude[noise_ind_start:]) # if system has less than 15 THz assume no artefacts at high frequencies and use til end
        return signal.amplitude/noise_floor

    noise_floor = np.mean(signal.amplitude[noise_ind_start:noise_ind_end])
    return signal.amplitude/noise_floor

def calc_t_meas(sample, reference):
    '''returns t_meas'''
    A = sample.amplitude / reference.amplitude
    phi = sample.phase - reference.phase
    return A, phi

def calc_complex_n(sample, reference):
    '''returns complex sample n_hat = n-jk as two arrays (n and k)'''

    A, phi = calc_t_meas(sample, reference)

    w = 2*np.pi*sample.frequency
    d = sample.thickness
    n_ref = reference.n
    
    n_sample = 1 - c*phi/d/w
    
    t_co = 4*n_sample*n_ref/(n_sample+n_ref)**2
    k_sample = -c/d/w*np.log(A/t_co)
    
    return n_sample, k_sample

def calc_FP(sample, reference, p=np.inf):
    d = sample.thickness
    w = sample.frequency
    n_s_complex = sample.n - 1j*sample.k
    n_r_complex = reference.n - 1j*reference.k

    X = ((n_s_complex - n_r_complex) / (n_s_complex + n_r_complex) * np.exp(-1j*w*d/c*n_s_complex))**2

    if type(p) == int:
        FP = np.full((len(X),), 0+0j)
        for i in range(p+1):
            FP += X**i
    else:
        FP = 1/(1-X)
    return FP

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
    global system_name
    system_name = var.get()
    sample_thickness = float(thickness_entry.get())*10**-6
    reference_n = float(reference_n_entry.get())
    root.quit()

def calc_theory_k(sample, T=300):
    A = 4.044*10**4
    B_0 = 1.391*10**5
    v_0 = 233
    f = sample.frequency
    n = sample.n
    h = 6.62607004*10**-34
    k = 1.38064852

    return 1/(2*n)*(A/f + (2*n*B_0*np.exp(h*c*v_0/k/T)/(4*np.pi*c*T*(np.exp(h*c*v_0/k/T)-1)**2*v_0**2))*f)

def set_skip_sep(system_name):
    if system_name == 'teraview':
        skip = 3
        sep = ','
    elif system_name == 'toptica':
        skip = 1
        sep = ','
    elif system_name == 'menlo':
        skip = 7
        sep = '\t'
    else:
        print(f'{system_name} not found.')
    return skip, sep

if __name__ == '__main__':
    
    # get user inputs
    root = tk.Tk()
    var = tk.StringVar()
    var.set('menlo')

    sample_filename, reference_filename = get_filenames()
    tk.Radiobutton(root, text='Teraview', variable=var, value='teraview').grid(row=0, column=0)
    tk.Radiobutton(root, text='TOPTICA', variable=var, value='toptica').grid(row=0, column=1)
    tk.Radiobutton(root, text='Menlo', variable=var, value='menlo').grid(row=1, column=0)
    
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
    skip, sep = set_skip_sep(system_name)
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
    second_peak_ps_sample_corrected = input('Input alternative FP peak location (ps) or press Enter if location found was acceptable:  ')
    if second_peak_ps_sample_corrected != '':
        second_peak_ps_sample = float(second_peak_ps_sample_corrected)
    
    sample_window = create_window(df_sample ,second_peak_ps_sample)
    ref_window = create_window(df_sample, second_peak_ps_sample, offset=main_peak_delay_ps)
    
    # calculate optical parameters
    sample = create_signal(df_sample, window=sample_window)
    reference = create_signal(df_reference, window=ref_window)  
    reference.n = reference_n
    sample.thickness = sample_thickness

    # estimate thickness with known estimate of average refractive index. For ice n~=1.81105
    # auto_thickness = c*10**-6/1.81105*(second_peak_ps_sample - main_peak_ps_sample)/2
    # use_auto_thickness = input(f'Estimated thicknes of ice is {auto_thickness} um; accept [y]/n ?:  ')
    # if use_auto_thickness != 'n':
    #     sample.thickness = auto_thickness * 10**-6
    
    # calculate optical parameters
    sample.n, sample.k = calc_complex_n(sample, reference)
    sample.absorption = 2*(2*np.pi*sample.frequency)/c*sample.k/100
    sample.permittivity = (sample.n**2 - sample.k**2) - 1j*2*sample.n*sample.k
     
    # calculate power spectra
    max = reference.amplitude.max()
    sample.power = calc_power(sample, max)
    reference.power = calc_power(reference, max)
    
    # calculate maximum absorption
    DR = calc_DR(reference)
    a_max = 2*np.log(DR*4*sample.n/(sample.n+1)**2)/(sample.thickness*100)

    # plots for display
    fig, axs = plt.subplots(2, 2)
    freq_range = np.arange(0, 5.5, 0.5)
    axs[0,0].plot(sample.frequency*10**-12, sample.n)
    
    axs[0,0].set_ylim([0,3])
    axs[0,0].set_xlabel('Frequency (THz)')
    axs[0,0].set_ylabel('n')
    axs[0,0].set_xticks(freq_range)
    axs[0,0].set_xlim([0.2, 5])

    axs[1,0].plot(sample.frequency*10**-12, sample.k)
    axs[1,0].set_xticks(freq_range)
    axs[1,0].set_xlim([0.2, 5])
    axs[1,0].set_xlabel('Frequency (THz)')
    axs[1,0].set_ylabel('k')
    axs[1,0].set_ylim([-0.1,0.4])

    axs[0,1].plot(sample.frequency*10**-12, sample.absorption, label='Sample')
    axs[0,1].plot(sample.frequency*10**-12, a_max, label='Max detectable')
    axs[0,1].set_xticks(freq_range)
    axs[0,1].set_xlim([0.2, 5])
    axs[0,1].set_ylim(bottom=0)
    axs[0,1].set_xlabel('Frequency (THz)')
    axs[0,1].set_ylabel('Absorption coefficient ($cm^{-1}$)')
    axs[0,1].legend()

    axs[1,1].plot(sample.frequency*10**-12, sample.power, label='Sample')
    axs[1,1].plot(reference.frequency*10**-12, reference.power, label='Reference') 
    axs[1,1].set_xlim([0, 10])
    axs[1,1].set_xlabel('Frequency (THz)')
    axs[1,1].set_ylabel('Power (dB)')
    axs[1,1].legend()
    
    plt.tight_layout()
    plt.show()

    # save results
    export_array = np.stack((sample.frequency*10**-12, sample.n, sample.absorption, a_max, sample.power, reference.power)).T[1:]
    export_df = pd.DataFrame(export_array, columns=['THz','n','a', 'max detectable', 'sample_power', 'ref_power'])

    try:
        with filedialog.asksaveasfile(mode='w', defaultextension=".csv") as file:
            export_df.to_csv(file.name, index=False)

    except AttributeError:
        # if user cancels save, filedialog returns None rather than a file object, and the 'with' will raise an error
        print("The user cancelled save")