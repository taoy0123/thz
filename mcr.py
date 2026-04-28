'''
Perform Multivariate Curve Resolution of spectra.
Takes inputs of spectral values in columns in .csv format. The top row contains numerical conditions such as temperatures (e.g., 298K). The left most column contains wavenumbers in cm-1.
Upper frequency limit set to 400 cm-1 (adjust as required).
Outputs four components and their contributions (adjust as required).
Residual spectra should be small amplitude and noise-like with no consistent peaks, otherwise add components to remove structure.
Residual norm should be roughly flat. If increasing with concentration or non-monotonic then model missing physics at high ionic strength or something changing (protonation).
Heatmap: vertical stripes mean missing spectral feature at a specific frequency. Horizontal patterns mean certain samples poorly fit. Bands changing with sample index mean evolving physics not captured.
Mean residual should be zero everywhere. If not then consistent under/over fitting at certain frequencies.
Created by Yu Heng Tao 09 Apr 2026
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import re
import statistics as stats

from tkinter.filedialog import askopenfilename, askopenfilenames
from tkinter.messagebox import showinfo
from tkinter import Tk, filedialog
import tkinter as tk

from pymcr.mcr import McrAR
from pymcr.constraints import ConstraintNonneg, Constraint

from scipy.optimize import nnls
from scipy.signal import savgol_filter

class ConstraintFixComponent(Constraint):
    def __init__(self, fixed_spectrum, index):
        self.fixed_spectrum = fixed_spectrum
        self.index = index

    def transform(self, ST):
        ST[self.index, :] = self.fixed_spectrum
        return ST

class ConstraintClosure(Constraint):
    def __init__(self, total=1.0): # total sums to 1
        self.total = total

    def transform(self, A):
        row_sums = A.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1  # avoid division by zero
        return A / row_sums * self.total

class ConstraintSmooth(Constraint):
    def __init__(self, window_length=11, polyorder=3): # window length typically between 9-11, polyorder 2 to 3
        self.window_length = window_length
        self.polyorder = polyorder

    def transform(self, A):
        return savgol_filter(
            A,
            window_length=self.window_length,
            polyorder=self.polyorder,
            axis=1
        )

def read_data(display, start_wavenum, end_wavenum):
    '''reads data into dataframe'''
    showinfo(message=display) 
    filename = askopenfilename()
    df = pd.read_csv(filename, sep=',', engine='python')
    reducedSpectra = df.loc[(df['wavenumber']>start_wavenum)&(df['wavenumber']<end_wavenum)]
    wavenumbers = reducedSpectra.iloc[:,0]
    data = reducedSpectra.iloc[:,1:]
    conditions = df.columns.values[1:]
    return wavenumbers, np.array(data.T), conditions

def normalize(v):
    return v / np.max(np.abs(v))

# Step 1 — Prepare inputs
# X: your measured spectra (glycine + HCl)
# shape: (n_samples, n_freq)

def save_result(df, parameter, index=True):
    '''saves parameter as csv file'''
    showinfo(message=f'Select save location of {parameter}') 
    try:
        with filedialog.asksaveasfile(mode='w', defaultextension=".csv") as file:
            df.to_csv(file.name, index=index)

    except TypeError:
    # if user cancels save, filedialog returns None rather than a file object, and the 'with' will raise an error
        print(f"The user cancelled saving {parameter}")


# asks custom frequency range
start_wavenum = float(input('Starting wavenumber (cm-1):  '))
end_wavenum = float(input('Ending wavenumber (cm-1):  '))
print(f'Performing MCR from {start_wavenum} to {end_wavenum} cm-1')

# import spectra
freq, X, hcl_conc = read_data('Select mixture file', start_wavenum, end_wavenum)
# Known spectra
_, S_hcl, _ = read_data('Select acid or base file', start_wavenum, end_wavenum)
_, S_gly, _ = read_data('Select amino acid file', start_wavenum, end_wavenum)

fix_hcl = ConstraintFixComponent(S_hcl, index=0)
fix_gly = ConstraintFixComponent(S_gly, index=1)

# Step 2 — Build initial Sᵀ
n_freq = X.shape[1]

# # random small guesses for unknown components
# np.random.seed(0)
# S_unk1 = np.random.rand(n_freq)

# initialise unknown component 1 using the residuals of nnls
A = np.vstack([S_hcl, S_gly]).T
sample = X[2]
coeffs, _ = nnls(A, sample)
fit = A @ coeffs
residual = sample - fit
S_unk1 = residual

# initialise unknown component 2
S_unk2 = sample - np.mean([S_hcl, S_gly], axis=0)

ST_init = np.vstack([
    normalize(S_hcl),
    normalize(S_gly),
    normalize(S_unk1),
    normalize(S_unk2)
])  # shape: (4, n_freq)


# Step 4 — Run MCR
mcr = McrAR(
    c_constraints=[
        ConstraintNonneg(),
        # ConstraintClosure(total=1.0)   # only if justified
    ],
    st_constraints=[
        ConstraintNonneg(),
        # ConstraintSmooth(window_length=9, polyorder=2),  # only if justified
        fix_hcl, 
        fix_gly
    ],
    max_iter=200
)

mcr.fit(X, ST=ST_init)

C = mcr.C_      # (n_samples, n_components)
ST = mcr.ST_    # (n_components, n_freq)

# Plot
# Example labels (edit as needed)
component_labels = ['HCl', 'Glycine', 'Unknown 1', 'Unknown 2']

# Create figure
fig, axes = plt.subplots(2, ST_init.shape[0], figsize=(16, 8))

# --- Top row: Component spectra ---
for i in range(ST_init.shape[0]):
    ax = axes[0, i]
    ax.plot(freq, ST[i])
    ax.set_title(f'{component_labels[i]} Spectrum')
    ax.set_xlabel('Frequency')
    ax.set_ylabel('Absorption')

# --- Bottom row: Concentrations vs HCl ---
for i in range(ST_init.shape[0]):
    ax = axes[1, i]
    ax.plot(hcl_conc, C[:, i], marker='o')
    ax.set_title(f'{component_labels[i]} Contribution')
    ax.set_xlabel('HCl Concentration (mM)')
    ax.set_ylabel('Contribution')

plt.tight_layout()
plt.show()

# Residuals
R = X - C @ ST
R = np.array(R)

# Plot residuals
# should be small amplitude and noise-like with no consistent peaks
# otherwise add components to remove structure

for i in range(R.shape[0]):
    plt.plot(freq, R[i], label=f'Sample {i}')

plt.xlabel('Frequency')
plt.ylabel('Residual')
plt.title('Residual Spectra')
plt.legend(hcl_conc)
plt.show()

# Plot residual norm
# should be roughly flat
# if increasing with concentration or non-monotonic, then model missing physics at high ionic strength or something changing

res_norm = np.linalg.norm(R, axis=1)
plt.plot(hcl_conc, res_norm, marker='o')
plt.xlabel('HCl Concentration (mM)')
plt.ylabel('Residual Norm')
plt.title('Fit Error vs HCl Concentration')
plt.show()

# Plot heatmap
# vertical stripes mean missing spectral feature at a specific frequency
# horizontal patterns mean certain samples poorly fit
# Bands changing with sample index mean evolving physics not captured

plt.imshow(R, aspect='auto')
plt.colorbar(label='Residual')
plt.xlabel('Frequency Index')
plt.ylabel('Sample Index')
plt.title('Residual Heatmap')
plt.show()

# Mean residual
# should be zero everywhere
# if not then consistent under/over fitting at certain frequencies
mean_residual = np.mean(R, axis=0)
plt.plot(freq, mean_residual)
plt.xlabel('Frequency')
plt.ylabel('Mean Residual')
plt.title('Average Residual Spectrum')
plt.show()

# Save results
save = input('Save results y/[n]? :   ')
if save == 'y':
    # construct dataframes
    df_spectra = pd.DataFrame(ST.T, columns=[
    'HCl', 'Gly', 'Unknown1'
    ])
    df_spectra['Frequency'] = freq

    df_conc = pd.DataFrame(C, columns=[
    'HCl', 'Gly', 'Unknown1'
    ])
    df_conc['HCl_concentration'] = hcl_conc

    save_result(df_spectra, 'spectra', index=False)
    save_result(df_conc, 'concentrations', index=False)