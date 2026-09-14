'''
Calculates the (de)protonated amino acid spectrum
Inputs: a_eff of amino acid solutions, their PCs and scores from PCA
Also requires a_eff of acid/base at 1M, and of amino acids without acid/base at 1M
Also need to input PKa or PKb of the amino acid
Can specify frequency range (must be the same as PCA)
Created by Yu Heng Tao 16 Apr 2026
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
from tkinter.filedialog import askopenfilename
from tkinter.messagebox import showinfo
from tkinter import filedialog
from scipy.linalg import svd

# Define the fitting function
def score_fit(conc_ab, A, B, conc_aa, Kab):
    x = conc_ab + conc_aa + Kab
    beta = (x-np.sqrt(-4*conc_ab*conc_aa+(x)**2))/2
    score = A * conc_ab + B * beta
    return score

def read_data(message, start_x=0, end_x=None):
    '''reads data into dataframe and extract data from range of interest'''
    showinfo(message=message) 
    filename = askopenfilename()
    df = pd.read_csv(filename, sep=',', engine='python')
    # if no upper limit of x specified, use the maximum
    if end_x is None:
        end_x = df.iloc[:, 0].max()
    # extract data from range of interest
    df_sliced = df.loc[(df.iloc[:, 0]>=start_x)&(df.iloc[:, 0]<=end_x)]
    x_values = df_sliced.iloc[:,0].to_numpy()
    y_values = df_sliced.iloc[:,1:].to_numpy()
    conditions = df.columns.values[1:].astype(float)
    return x_values, y_values.T, conditions

def save_result(df, message, index=True):
    '''saves parameter as csv file'''
    showinfo(message=message) 
    try:
        with filedialog.asksaveasfile(mode='w', defaultextension=".csv") as file:
            df.to_csv(file.name, index=index)

    except TypeError:
    # if user cancels save, filedialog returns None rather than a file object, and the 'with' will raise an error
        print("The user cancelled save")

# specify concentration of the amino acid (M)
conc_aa = 1

# specify PKa (carboxyl) or PKb (amino) of the amino acid
# Glycine: PKa = 2.34, PKb = 9.6
# Serine: PKa = 2.21, PKb = 9.15
# Proline: PKa = 1.99, PKb = 10.96
PKab = 10.96
Kab = 10 ** (-PKab)

# specify frequency range (must be the same as PCA)

try:
    start_wavenum = float(input('Starting wavenumber (cm-1):  '))
# if user input invalid set to beginning of input data
except ValueError:
    start_wavenum = 0
try:
    end_wavenum = float(input('Ending wavenumber (cm-1):  '))
# if user input invalid set 400 cm-1 as upper limit
except ValueError:
    end_wavenum = 400
print(f'Performing SVD from {start_wavenum} to {end_wavenum} cm-1')

wavenumbers, data, conditions = read_data('Select data file', start_x=start_wavenum, end_x=end_wavenum)
# dataMatrix = np.array(data)
U,s,Vh = svd(data)

population = U
weight = s/s.sum()

sigma = np.zeros(data.shape)
for i in range(min(sigma.shape)):
    sigma[i,i] = s[i]

spectrum = np.dot(sigma, Vh)

basis_spectra = []
populations = []

# plot outputs
rank = 4
fig, axes = plt.subplots(nrows=2, ncols=rank)
for i in range(rank):
    # plot basis spectra
    # flip spectrum and score if negative
    if np.sum(population[:,i]<0) >= np.sum(population[:,i]>0):
        spectrum[i,:] = -spectrum[i,:]
        population[:,i] = -population[:,i]
    axes[0,i].plot(wavenumbers,spectrum[i,:])
    axes[0,i].set_title(f'Principal Component {i+1} ({weight[i]*100:.2f}%)')
    axes[0,i].set(xlabel='Wavenumber (cm$^{-1}$)', ylabel='Signal intensity (a.u.)')
    basis_spectra.append(spectrum[i,:])
    # plot scores
    axes[1,i].plot(conditions, population[:,i], marker='o', linestyle='-')
    axes[1,i].set(xlabel='Conditions', ylabel='Population (a.u.)')
    populations.append(population[:,i])
plt.show()


PCs = np.array(basis_spectra)
PC1 = PCs[0]
PC2 = PCs[1]

scores = np.array(populations)
score_1 = scores[0]
score_2 = scores[1]

# import scores and acid/base concentrations
conc_ab = conditions / 1000
# conc_ab = conditions

# Fit score 1 with fitting function
parameters_1, covariance_1 = curve_fit(
    lambda conc_ab, A, B: score_fit(conc_ab, A, B, conc_aa, Kab),
    conc_ab, 
    score_1
)

# Fit score 2 with fitting function
parameters_2, covariance_2 = curve_fit(
    lambda conc_ab, A, B: score_fit(conc_ab, A, B, conc_aa, Kab),
    conc_ab, 
    score_2
)

A = [parameters_1[0], parameters_2[0]]
B = [parameters_1[1], parameters_2[1]]

fit_y1 = score_fit(conc_ab, parameters_1[0], parameters_1[1], conc_aa, Kab)
fit_y2 = score_fit(conc_ab, parameters_2[0], parameters_2[1], conc_aa, Kab)

fig, axes = plt.subplots(nrows=1, ncols=2)

axes[0].set_title('Score of PC1')
axes[0].plot(conc_ab, score_1, 'o', label='data')
axes[0].plot(conc_ab, fit_y1, '-', label='fit')
axes[0].set(xlabel='Acid/Base concentration (M)', ylabel='Score')

axes[1].set_title('Score of PC2')
axes[1].plot(conc_ab, score_2, 'o', label='data')
axes[1].plot(conc_ab, fit_y2, '-', label='fit')
axes[1].set(xlabel='Acid/Base concentration (M)', ylabel='Score')

plt.legend()
plt.show()

# Save results
save_result(pd.DataFrame(
    np.column_stack((wavenumbers, PCs.T))
    ), 'Select save location of PCs', index=False)
save_result(pd.DataFrame(
    np.column_stack((conc_ab, np.array(populations).T))
    ), 'Select save location of scores', index=False)
save_result(pd.DataFrame(
    np.column_stack((conc_ab, fit_y1.T, fit_y2.T))
    ), 'Select save location of fits', index=False)
save_result(pd.DataFrame(np.array(weight).T), 'Select save location of weights', index=False)

# # reconstruct using score 1 and 2
# score_1 = score_1[:, None]
# score_2 = score_2[:, None]
# reconstructed_spectra = PC1 * score_1 + PC2 * score_2

# # reconstruct using fitted scores
# fit_y1 = fit_y1[:, None]
# fit_y2 = fit_y2[:, None]
# reconstructed_spectra = PC1 * fit_y1 + PC2 * fit_y2

# reconstruct spectra  using A and B
conc_ab = conc_ab[:, None]
x = conc_ab + conc_aa + Kab
beta = (x-np.sqrt(-4*conc_ab*conc_aa+(x)**2))/2
reconstructed_spectra = (A[0] * PC1 + A[1] * PC2) * conc_ab + (B[0] * PC1 + B[1] * PC2) * beta
reconstructed_spectra_ab = (A[0] * PC1 + A[1] * PC2)

# import measured acid/base spectrum
_, a_ab, _ = read_data('Select acid/base spectrum file', start_x=start_wavenum, end_x=end_wavenum)
a_ab = a_ab[0]

# check reconstructed spectra vs original
plt.title('Effective absorption: data vs reconstructed')
plt.plot(wavenumbers, data[2], 'o', label='Measured (1M)')
plt.plot(wavenumbers, reconstructed_spectra[2], '-', label='Reconstructed (1M)')
plt.plot(wavenumbers, a_ab,  label='Measured HCl or NaOH (1M)')
plt.plot(wavenumbers, reconstructed_spectra_ab, label='Reconstructed HCl or NaOH (1M)')
plt.xlabel('Wavenumber (cm-1)')
plt.ylabel('Effective absorption coefficient')
plt.legend()
plt.show()

# save reconstructed acid or base spectrum
df_da_ab = pd.DataFrame({
    "wavenumber": wavenumbers,
    "HCl/NaOH": reconstructed_spectra_ab
})
save_result(df_da_ab, 'Select save location of HCl/NaOH spectrum', index=False)


# calculate protonated/deprotonated amino acid spectrum
# import measured amino acid spectrum
_, a_aa, _ = read_data('Select amino acid spectrum file', start_x=start_wavenum, end_x=end_wavenum)
a_aa = a_aa[0]


# calculate difference spectra
da = (B[0] * PC1 + B[1] * PC2)

# #calculate HCl or NaOH contributions
# da_ab = (A[0]*PC1 + A[1]*PC2)

# calculate protonated/deprotonated amino acid spectrum
a_charged_aa = da + a_aa + a_ab

# plot protonated/deprotonated amino acid spectrum
plt.title('Effective absorption: zwitterion vs (de)protonated')
plt.plot(wavenumbers, a_aa, label='zwitterion')
plt.plot(wavenumbers, a_charged_aa, label='(de)protonated')
plt.xlabel('Wavenumber (cm-1)')
plt.ylabel('Effective absorption')
plt.legend()
plt.show()

# Save results
# save (de)protonated spectrum
df_charged_aa = pd.DataFrame({
    "wavenumber": wavenumbers,
    "(de)protonated": a_charged_aa
})
save_result(df_charged_aa, 'Select save location of (de)protonated spectrum', index=False)

