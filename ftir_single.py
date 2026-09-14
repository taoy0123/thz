'''
Computes absorbance difference between two samples measured with Bruker Vertex 80v in transmission.
Inputs single beam intensities recorded as OPUS files and exported as text by OPUS file reader.
Created by Yu Heng Tao 18 Nov 2022
Acknowledgement to Simon Schulke and Gerhard Schwaab.
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

from tkinter.filedialog import askopenfilename
from tkinter.messagebox import showinfo
from tkinter import Tk, filedialog

from scipy.signal import find_peaks, butter, freqs, savgol_filter, sosfiltfilt, argrelmin
from scipy.fft import rfft, rfftfreq, irfft
from scipy import interpolate, fftpack

def read_data(filename, sep='\t'):
    df = pd.read_csv(filename, sep=sep, names=['wavenumber', 'intensity'], engine='python')
    df = df.sort_values(by=['wavenumber'])
    return df

def get_filename():
    '''get names of input files'''
    showinfo(message='Select sample file')  
    sample_filename = askopenfilename()
    # switch out the below lines for pre-defined reference file
    showinfo(message='Select reference file')
    reference_filename = askopenfilename()
    # reference_filename = 'C:/Users/Tao/Documents/lab_results/2026-04-23_Hank/04_Water_13um_THz_RT_MM-1p5mm_Globar_2cm-1_128Scans_ev-v.0.txt'

    return sample_filename, reference_filename


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
    # if there is no peak in the range, don't apply any window
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
    # 150 cm-1 peak height in vapor spectrum is 30.42173
    scaling_factor = vapor_amount / 30.42173
    residual_vapor = scaling_factor * vapor_spectrum[:len(xf)]

    # subtract the residual vapor spectrum from the signal
    if xf[vapor_max_index] < xf[vapor_zero_index]:
        corrected_spectrum = xf + residual_vapor
    else:
        corrected_spectrum = xf - residual_vapor

    # # plot results
    # plt.plot(f, xf)
    # plt.plot(f, corrected_spectrum)
    # plt.plot(f, residual_vapor)
    # ax = plt.gca()
    # ax.set_xlim([50,600])
    # ax.set_xlabel('Wavenumbers (cm$^{-1}$)')
    # ax.set_ylabel(r'$\Delta$ A')
    # plt.legend(['Raw', 'Vapor corrected', 'Residual vapor'])
    # plt.show()

    return residual_vapor, corrected_spectrum

def lowpass(data: np.ndarray, cutoff: float, sample_rate: float, poles: int = 5):
    sos = butter(poles, cutoff, 'lowpass', fs=sample_rate, output='sos')
    filtered_data = sosfiltfilt(sos, data)
    return filtered_data

def cell_thickness():
    # determine cell thickness from etalons in the transmitted spectrum of empty cell measurement

    # load empty cell intensities
    showinfo(message='Select empty cell file')  
    empty_cell_filename = askopenfilename()
    empty_cell_df = pd.read_csv(empty_cell_filename, sep='\t', names=['wavenumber', 'intensity'], engine='python')
    empty_cell_df = empty_cell_df.sort_values(by=['wavenumber'])

    # load background intensities
    showinfo(message='Select background file') 
    background_filename = askopenfilename()
    background_df = pd.read_csv(background_filename, sep='\t', names=['wavenumber', 'intensity'], engine='python')
    background_df = background_df.sort_values(by=['wavenumber'])

    # slice background intensities to the same frequency range as empty cell's
    starting_wavenumber = empty_cell_df['wavenumber'].iloc[0]
    ending_wavenumber = empty_cell_df['wavenumber'].iloc[-1]
    background_df = background_df[background_df['wavenumber'].between(starting_wavenumber, ending_wavenumber)]

    f = empty_cell_df['wavenumber'].values
    xf = empty_cell_df['intensity'].values

    bf = background_df['wavenumber'].values
    bxf = background_df['intensity'].values

    # THz = f*0.03

    # compute difference absorbance
    delta_A = -np.log(bxf/xf)
    # replace NaN values with interpolated values
    delta_A = pd.Series(delta_A).interpolate().to_numpy()

    # convert nan to zero. nan due to -ve intensity
    np.nan_to_num(delta_A, copy=False)

    # # Remove water vapor lines in empty cell and background measurement
    # # Comment out if using mid IR data (vapor spectrum in THz only)
    # residual_vapor, delta_A_vapor_corrected = vapor_correction(f, delta_A, vapor_intensity)

    # remove etalons caused by the cell windows (MIR)
    delta_A_smoothed = lowpass(delta_A, 0.05, 1)
    delta_A_smoothed = savgol_filter(delta_A_smoothed, 150, 2)

    # locate minima in delta A spectrum
    minima_array = argrelmin(delta_A_smoothed, order=10)[0]

    
    # look up corresponding wavenumbers of minima
    minima_wavenumbers = [f[i] for i in minima_array]
    # only include minimas within 4000-7000 cm-1, where the etalons are the most apparent. If measured in THz, use the entire frequency range (i.e., comment out the below line)
    # minima_wavenumbers = [minima_wavenumber for minima_wavenumber in minima_wavenumbers if minima_wavenumber > 4000 and minima_wavenumber < 6800]

    # plot
    plt.plot(f, delta_A)
    plt.plot(f, delta_A_smoothed)
    plt.vlines(minima_wavenumbers, min(delta_A_smoothed), max(delta_A_smoothed), colors="red", linestyles="dashed")
    plt.xlabel("Wavenumber (cm-1)")
    plt.ylabel("Intensity (arb. u.)")
    plt.legend(['raw', 'smoothed'])
    plt.show()

    # calculate the average distance between minimas
    average_period = (max(minima_wavenumbers) - min(minima_wavenumbers)) / (len(minima_wavenumbers) - 1)

    # The reciprocal of the period is twice the cell thickness. calculate cell thickness
    cell_thickness_um_by_minima = (1/average_period)*1e4 / 2 

    print(f'The empty cell thickness by minima is {cell_thickness_um_by_minima} \u03bcm.')

    # fourier transform
    dF = f[2] - f[1]
    xt = rfft(delta_A_smoothed)
    t = rfftfreq(len(delta_A_smoothed), dF)

    # exclude beginning of fourier transform in peak search (<0.0015)
    start_index = np.where(t>10*2/1e4)[0][0]
    amp_xt = abs(xt)
    max_amp_xt = amp_xt[start_index:].max()
    peak_index = np.where(amp_xt==max_amp_xt)[0][0]

    cell_thickness_um_by_fourier = (t[peak_index]*1e4)/2
    cell_thickness_um_resolution = ((t[1] - t[0])*1e4)/2
    print(f'The empty cell thickness by Fourier Transform is {cell_thickness_um_by_fourier} \u03bcm')

    # choose which method to use
    # cell_thickness_um = cell_thickness_um_by_minima
    cell_thickness_um = cell_thickness_um_by_fourier

    min_cell_thickness_um = cell_thickness_um - (cell_thickness_um_resolution/2)
    max_cell_thickness_um = cell_thickness_um + (cell_thickness_um_resolution/2)
    print(f'The empty cell thickness is between {round(min_cell_thickness_um, 2)} and {round(max_cell_thickness_um, 2)} \u03bcm')
    
    # # plot Fourier Transform of intensity
    # plt.bar(t, amp_xt, width=0.0001, tick_label=np.round((t*1e4)/2,1))
    # plt.xlabel("Cell thickness (\u03bcm)")
    # plt.ylabel("Intensity (arb. u.)")
    # plt.xlim(5*2/1e4, 40*2/1e4)
    # plt.ylim(0, max_amp_xt*1.1)
    # plt.xticks(fontsize=7, rotation=90)
    # plt.legend(['smoothed'])
    # plt.show()

    return min_cell_thickness_um, max_cell_thickness_um

def asymmetric_least_squares(y, x, lam=1e6, p=0.01, niter=100,
                             xmin=50, xmax=600):
    """
    Asymmetric least squares baseline fitting, restricted to a wavenumber range.

    Parameters
    ----------
    y : array
        Spectrum intensities (1D).
    x : array
        Wavenumber axis (1D, same length as y).
    lam : float
        Smoothness parameter. Bigger the lam, flatter the baseline.
    p : float
        Asymmetry parameter (0.001–0.01). Smaller the p, lower the baseline.
    niter : int
        Number of iterations.
    xmin, xmax : float
        Wavenumber range to apply baseline fitting.

    Returns
    -------
    baseline : array
        Full-length baseline (same size as y). Baseline fitted only in [xmin, xmax].
    """

    # mask for the fitting region
    mask = (x >= xmin) & (x <= xmax)
    y_fit = y[mask]
    L = len(y_fit)

    # Build second difference matrix
    D = np.zeros((L-2, L))
    for i in range(L-2):
        D[i, i]   = 1
        D[i, i+1] = -2
        D[i, i+2] = 1
    DTD = D.T @ D

    w = np.ones(L)
    for _ in range(niter):
        W = np.diag(w)
        Z_fit = np.linalg.solve(W + lam * DTD, w * y_fit)
        w = p * (y_fit > Z_fit) + (1 - p) * (y_fit < Z_fit)

    # Build full baseline
    baseline = np.zeros_like(y)
    baseline[mask] = Z_fit
    return x[mask], baseline[mask]

def remove_outliers(Delta_A, frequency):
    Delta_A = pd.Series(Delta_A, index=frequency)
    mask = (Delta_A.index >= 50) & (Delta_A.index <= 600)
    # --- Detect outliers ---
    window_size = 11
    rolling_median = Delta_A.rolling(window_size, center=True).median()
    residuals = Delta_A - rolling_median
    threshold = 15 * np.median(np.abs(residuals[mask].dropna()))
    outliers = (residuals.abs() > threshold) & mask
    outlier_idx = Delta_A.index[outliers]

    # remove outliers
    new_Delta_A = Delta_A.copy()
    for idx in Delta_A.index[outliers]:
        pos = Delta_A.index.get_loc(idx)
        if pos > 0 and pos < len(Delta_A) - 1:
            left = Delta_A.iloc[pos - 3]
            right = Delta_A.iloc[pos + 3]
            new_Delta_A.iloc[pos] = (left + right) / 2
    
    # # --- Plot without outliers ---
    # plt.plot(frequency, Delta_A.values, label="Data", color="green")
    # plt.plot(frequency, new_Delta_A.values, label="Data", color="blue")
    # plt.scatter(outlier_idx, Delta_A[outliers], color="red", label="Outliers")
    # for xi, yi in zip(outlier_idx, Delta_A[outliers]):
    #     plt.text(xi, yi, f"{xi:.1f}", fontsize=8, color="red", ha="center", va="bottom")
    # ax = plt.gca()
    # ax.set_xlim([50,600])
    # plt.legend(['original','without outliers', 'outliers'])
    # plt.xlabel("Wavenumber")
    # plt.ylabel("Delta A")
    # plt.show()

    return new_Delta_A

def remove_fringes_halfperiod(wavenumbers, absorbance, period=None, plot_fft=False):
    """
    Remove etalon fringes from an FTIR absorption spectrum by shifting
    the spectrum by half a fringe period and averaging.

    If period is not provided, it is estimated automatically from the
    200–300 cm⁻¹ spectral region, assuming the fringes are most visible there.

    Parameters
    ----------
    wavenumbers : array-like
        1D array of wavenumber values (monotonic, either increasing or decreasing).
    absorbance : array-like
        1D array of absorbance (or transmission) values corresponding to wavenumbers.
    period : float, optional
        Known fringe period in cm⁻¹. If None, it will be estimated automatically.
    plot_fft : bool, optional
        If True, plots the FFT amplitude spectrum used to determine the period.

    Returns
    -------
    wn_trunc : ndarray
        Truncated wavenumber axis (excluding edges affected by the shift).
    abs_no_fringes : ndarray
        Absorbance values with fringes removed.
    detected_period : float
        Estimated or provided fringe period (in cm⁻¹).
    """

    wavenumbers = np.asarray(wavenumbers)
    absorbance = np.asarray(absorbance)

    # Ensure increasing wavenumber order
    if wavenumbers[0] > wavenumbers[-1]:
        wn = wavenumbers[::-1]
        ab = absorbance[::-1]
        reversed_output = True
    else:
        wn = wavenumbers
        ab = absorbance
        reversed_output = False

    # --- 1. Determine fringe period if not provided ---
    if period is None:
        # Focus only on the custom 200-300 cm⁻¹ range
        mask_window = (wn >= 520) & (wn <= 540)
        if np.sum(mask_window) < 10:
            raise ValueError("Not enough data points in cm⁻¹ region for FFT analysis.")

        wn_sub = wn[mask_window]
        ab_sub = ab[mask_window]

        # Detrend baseline in that region
        ab_detrended = ab_sub - np.polyval(np.polyfit(wn_sub, ab_sub, 2), wn_sub)
        y = ab_detrended - np.mean(ab_detrended)

        # Compute FFT
        dnu = np.mean(np.diff(wn_sub))
        freqs = fftpack.rfftfreq(len(wn_sub), d=dnu)  # cycles per cm⁻¹
        spectrum = np.abs(fftpack.rfft(y))

        # Expected fringe period ~4 cm⁻¹ => frequency ~ 1/4 = 0.25 cycles/cm⁻¹
        # Search reasonable band (2–10 cm⁻¹ → 0.1–0.5 cycles/cm⁻¹)
        f_min, f_max = 1/10.0, 1/2.0
        band_mask = (freqs >= f_min) & (freqs <= f_max)

        if not np.any(band_mask):
            raise ValueError("Frequency grid too coarse for 2–10 cm⁻¹ period search range.")

        idx = np.argmax(spectrum[band_mask])
        fringe_freq = freqs[band_mask][idx]
        period = 1.0 / fringe_freq

        if plot_fft:
            import matplotlib.pyplot as plt
            plt.plot(freqs, spectrum, label="FFT amplitude")
            plt.axvspan(f_min, f_max, color='orange', alpha=0.2, label="Search band (2–10 cm⁻¹)")
            plt.xlabel("Frequency (cycles per cm⁻¹)")
            plt.ylabel("Amplitude (a.u.)")
            plt.title(f"Detected fringe period ≈ {period:.2f} cm⁻¹ (from 200–300 cm⁻¹ region)")
            plt.legend()
            plt.show()

    # --- 2. Interpolate and shift by half-period ---
    interp = interpolate.interp1d(
        wn, ab, kind='linear', bounds_error=False, fill_value=np.nan
    )
    ab_shifted = interp(wn + period / 2.0)

    # --- 3. Compute average where both data exist ---
    mask = ~np.isnan(ab_shifted)
    wn_trunc = wn[mask]
    abs_no_fringes = 0.5 * (ab[mask] + ab_shifted[mask])

    # --- 4. Reverse back if original was decreasing ---
    if reversed_output:
        wn_trunc = wn_trunc[::-1]
        abs_no_fringes = abs_no_fringes[::-1]

    return wn_trunc, abs_no_fringes, period

if __name__ == '__main__':

    # read in water vapor data (valid up to 400 cm-1)
    target_path = os.path.join(os.path.dirname(__file__), 'WaterVapor_400.csv')
    vapor_df = read_data(target_path, sep=',')   
    vapor_frequency = vapor_df['wavenumber'].values
    vapor_intensity = vapor_df['intensity'].values

    # # calculate thickness
    # min_sample_thickness_um, max_sample_thickness_um = cell_thickness()
    # sample_thickness_um = (min_sample_thickness_um + max_sample_thickness_um) / 2

    ## 0528 Con A (3mg/ml)
    ## 0514 water
    ## 0709 BSA (1mg/ml)
    ## 0716 water
    # sample_thickness_um = 10.45667558127589
    
    # 0425, 0509 Con A (2mg/ml)
    # 0521 Con A (3mg/ml)
    # 0602 water
    # 0626 BSA (3mg/ml)
    # sample_thickness_um = 11.110217805105632

    # # 0611 Con A (3mg/ml)
    # sample_thickness_um = 11.110896492114588

    # # 0704 BSA (2mg/ml)
    # sample_thickness_um = 11.763760028935376

    # # 0804 Proline (0.1M)
    # # 0818 Proline (1M)
    # sample_thickness_um = 15.386361837798383

    # # 0829 Glycine
    # sample_thickness_um = 15.686104320492694

    # # 0915 Proline
    # # 0916 BSA 50mg/ml
    # # 0930 Serine
    # # 1201 Water
    # sample_thickness_um = 16.338555392419778

    # 1014 Cysteine 0.65M
    # sample_thickness_um = 19.606266470903734

    # 1120 BSA 1mg
    # sample_thickness_um = 17.64563982381336
    
    # 0115 ConA 50mg
    # sample_thickness_um = 18.299182039510153

    # 0318 HCl NaOH 0.3-2.5M
    # sample_thickness_um = 15.685013178349386

    # 0423 NaOH 0.3-2.5M
    # sample_thickness_um = 16.338555394113946

    # # 0719 BSA and water
    # sample_thickness_um = 20.913350904465847

    # # 0722 BSA and water
    # # 0730 Water and BSA
    # sample_thickness_um = 20.25980868870129

    # # 0806 Water and BSA 15mg/ml
    # sample_thickness_um = 19.60626647293673

    # 0807 Water and BSA 20mg/ml water thickness
    # BSA thickness 22.10174761105815
    # sample_thickness_um = 20.25980868870129

    # ATR
    sample_thickness_um = 10000


    sample_thickness_cm = sample_thickness_um * 1e-4 
    # min_sample_thickness_cm = min_sample_thickness_um * 1e-4 
    # max_sample_thickness_cm = max_sample_thickness_um * 1e-4 
    # sample_thickness_cm = 1 

    # read in measured data
    sample_filename, reference_filename = get_filename()
    sample_df = read_data(sample_filename)
    reference_df = read_data(reference_filename)

    frequency = sample_df['wavenumber'].values
    sample_intensity = sample_df['intensity'].values
    reference_intensity = reference_df['intensity'].values

    # compute difference absorbance
    delta_A = -np.log(sample_intensity/reference_intensity)
    # replace NaN values with interpolated values
    delta_A = pd.Series(delta_A).interpolate().to_numpy()


    # filter water vapor
    # residual_vapor, delta_A_corrected = vapor_correction(frequency, delta_A, vapor_intensity)
    delta_A_corrected = delta_A

    # # detect and remove outliers
    # delta_A_without_corrected = remove_outliers(delta_A_corrected, frequency)

    # remove window fringes
    filtered_frequency, delta_A_smoothed, detected_period = remove_fringes_halfperiod(frequency, delta_A_corrected, period=3.86, plot_fft=True)

    # # smooth spectrum using AsLS
    # # Bigger lam, flatter the baseline. Smaller p, lower the baseline. Asymmetry parameter (0.001–0.01)
    # filtered_frequency, delta_A_smoothed = asymmetric_least_squares(delta_A_corrected, frequency, lam=3e4, p=0.3)

    # # smooth spectrum using Savitzky-Golay filter, window size 80 and polynomial order 3 by trial and error
    # delta_A_smoothed = savgol_filter(delta_A_corrected, 80, 3)  

    # # plot difference absorbance and check the effect of smoothing
    # mask = (frequency >= 50) & (frequency <= 600)
    # plt.plot(frequency[mask], delta_A_corrected[mask])
    # plt.plot(filtered_frequency, delta_A_smoothed)
    # ax = plt.gca()
    # ax.set_xlim([50,600])
    # # ax.set_ylim([-0.5,0.5])
    # ax.set_xlabel('Wavenumbers (cm$^{-1}$)')
    # ax.set_ylabel(r'$\Delta$ A')
    # plt.legend(['Vapour Corrected', 'Fringes removed'])
    # plt.show()

     
    # # if not smoothing
    # delta_A_smoothed = delta_A_corrected

    
    delta_a_co = delta_A_smoothed / sample_thickness_cm
    # min_delta_a_co = delta_A_smoothed / max_sample_thickness_cm
    # max_delta_a_co = delta_A_smoothed / min_sample_thickness_cm

    # plus_error = max_delta_a_co - delta_a_co
    # minus_error = delta_a_co - min_delta_a_co

    # save results
    export_array = np.stack((filtered_frequency, delta_a_co)).T[1:]
    export_df = pd.DataFrame(export_array, columns=['wavenumber', 'delta a_co'])

    try:
        # # for removing BG in ATR
        # with filedialog.asksaveasfile(mode='w', defaultextension=".txt") as file:
        #     export_df.to_csv(
        #         file.name,
        #         sep='\t',
        #         index=False,
        #         header=False
        #     )
        with filedialog.asksaveasfile(mode='w', defaultextension=".csv") as file:
            export_df.to_csv(file.name, index=False)

    except TypeError:
        # if user cancels save, filedialog returns None rather than a file object, and the 'with' will raise an error
        print("The user cancelled save")