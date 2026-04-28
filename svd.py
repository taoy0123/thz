'''
Perform singular value decomposition of spectra
Takes inputs of spectral values in columns in .csv format. The top row contains numerical conditions such as temperatures (e.g., 298K). The left most column contains wavenumbers in cm-1
Upper frequency limit set to 400 cm-1 (adjust as required)
Outputs the four most significant components and their scores
Presents outputs as basis spectra, populations, and weights. Populations are dot product of singular values and right singular vectors. Weights are singular values in percentages.
Created by Yu Heng Tao 05 Jul 2023
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import re
import statistics as stats

from tkinter.filedialog import askopenfilename, askopenfilenames
from tkinter.messagebox import showinfo
from tkinter import filedialog

from scipy.linalg import svd

def read_data():
    '''reads data into dataframe'''
    showinfo(message='Select data file') 
    filename = askopenfilename()
    df = pd.read_csv(filename, sep=',', engine='python')
    return df

def save_result(df, parameter, index=True):
    '''saves parameter as csv file'''
    showinfo(message=f'Select save location of {parameter}') 
    try:
        with filedialog.asksaveasfile(mode='w', defaultextension=".csv") as file:
            df.to_csv(file.name, index=index)

    except TypeError:
    # if user cancels save, filedialog returns None rather than a file object, and the 'with' will raise an error
        print(f"The user cancelled saving {parameter}")

df = read_data()

# asks custom frequency range
try:
    start_wavenum = float(input('Starting wavenumber (cm-1):  '))
# if user input invalid set to beginning of input data
except ValueError:
    start_wavenum = min(df['wavenumber'])
try:
    end_wavenum = float(input('Ending wavenumber (cm-1):  '))
# if user input invalid set 400 cm-1 as upper limit
except ValueError:
    end_wavenum = 400
print(f'Performing SVD from {start_wavenum} to {end_wavenum} cm-1')


reducedSpectra = df.loc[(df['wavenumber']>=start_wavenum)&(df['wavenumber']<=end_wavenum)]
# store spectra and wavenumbers separately
wavenumbers = reducedSpectra.iloc[:,0]
data = reducedSpectra.iloc[:,1:]

# extract conditions (e.g., temperatures) from the header of each column
# skip first column (i.e., wavenumbers)
conditions = df.columns.values[1:]
# remove non-digits from each condition
for i,condition in enumerate(conditions):
    conditions[i] = float(re.sub(r'[^\d.]', '', condition))

# perform svd
dataMatrix = np.array(data.T)
U,s,Vh = svd(dataMatrix)

# # reconstruct sigma
# # set output matrix dimensions
# m = len(U)
# n = len(Vh)
# # construct the diagonal matrix from the singular values
# # create a matrix filled with zeroes
# sigma = np.zeros((m, n))
# # insert singular values
# for i in range(min(m, n)):
#     # reconstruct
#     sigma[i, i] = s[i]


# # Reconstruct spectra
# reconstructed_spectra = pd.DataFrame(np.dot(U, np.dot(sigma, Vh)))
# reconstructed_spectra.columns = conditions
# reconstructed_spectra.index = wavenumbers
# reconstructed_spectra.index.name = "wavenumbers"
# data.index = wavenumbers
# data.index.name = "wavenumbers"
# # # Subtract reconstructed spectra from original and calculate absolute error
# errors = abs(data - reconstructed_spectra)

# print(f'The maximum error percentage is {errors.to_numpy().max()}')

# save_result(reconstructed_spectra, 'Reconstructed spectra')
# save_result(errors, 'Errors')

population = U
# store singular values
weight = s/s.sum()
# convert singular values to matrix for dot product calculation
sigma = np.zeros(dataMatrix.shape)
for i in range(min(sigma.shape)):
    sigma[i,i] = s[i]
spectrum = np.dot(sigma, Vh)

# initialise lists for export
spectra = [wavenumbers]
populations = [conditions]

# plot outputs
rank = 4
fig, axes = plt.subplots(nrows=2, ncols=rank)
for i in range(rank):
    # plot basis spectra
    axes[0,i].plot(wavenumbers,spectrum[i,:])
    axes[0,i].set_title(f'Principal Component {i+1} ({weight[i]*100:.2f}%)')
    axes[0,i].set(xlabel='Wavenumber (cm$^{-1}$)', ylabel='Signal intensity (a.u.)')
    spectra.append(spectrum[i,:])
    # plot scores
    axes[1,i].plot(conditions, population[:,i], marker='o', linestyle='-')
    axes[1,i].set(xlabel='Temperatures ($^{o}$C)', ylabel='Population (a.u.)')
    populations.append(population[:,i])
plt.show()

# save results
# construct dataframes
export_spectra = pd.DataFrame(np.array(spectra).T)
export_populations = pd.DataFrame(np.array(populations).T)
# export_weights = pd.DataFrame(weight.T, columns=['weight normalised'])
# export_sigmas = pd.DataFrame(s.T, columns=['sigma'])

save_result(export_spectra, 'PCs', index=False)
save_result(export_populations, 'scores', index=False)
# save_result(export_weights, 'weights', index=False)
# save_result(export_sigmas, 'sigmas', index=False)