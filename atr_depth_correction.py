from tkinter.filedialog import askopenfilename
from tkinter.messagebox import showinfo
from tkinter import Tk, filedialog

import numpy as np
from scipy.interpolate import interp1d
import pandas as pd
import matplotlib.pyplot as plt
import os

def get_filename(prompt_message):
    '''get names of input files'''
    showinfo(message=prompt_message)  
    target_filename = askopenfilename()

    return target_filename

def read_water_data(filename, sep='\t'):
    df = pd.read_csv(filename, sep=sep, header=0, names=['wavenumber', 'n', 'k'], engine='python')
    df = df.sort_values(by=['wavenumber'])
    return df

def get_water_nk(wavenumbers, df_water_n_k):
    """
    Interpolate water n and k onto arbitrary wavenumbers.
    """

    wn_data = df_water_n_k['wavenumber'].values
    n_data = df_water_n_k['n'].values
    k_data = df_water_n_k['k'].values


    wavenumbers = np.asarray(wavenumbers, dtype=float)

    n_interp = interp1d(
        wn_data,
        n_data,
        kind="linear",
        bounds_error=True
    )

    k_interp = interp1d(
        wn_data,
        k_data,
        kind="linear",
        bounds_error=True
    )

    n = n_interp(wavenumbers)
    k = k_interp(wavenumbers)

    return n, k

def complex_qz(wavenumbers, n, k,
               n_crystal=2.4,
               theta_deg=45):

    theta = np.deg2rad(theta_deg)

    # Complex refractive index
    n_complex = n + 1j * k

    # Complex normal wavevector
    qz = (
        2 * np.pi * wavenumbers *
        np.sqrt(
            n_complex**2 -
            (n_crystal * np.sin(theta))**2
        )
    )

    return qz

def effective_depth_from_qz(qz):

    gamma = np.imag(qz)

    d_eff_cm = 1 / (2 * gamma)

    return d_eff_cm

if __name__ == '__main__':

  # Load water refractive index 
  water_n_k_filename = get_filename('Select water refractive index file')
  df_water_n_k = read_water_data(water_n_k_filename)

  df = pd.read_csv(get_filename('Select target wavenumbers file'), header=0, names=['wavenumber'], usecols=[0], engine='python')
  df = df.sort_values(by=['wavenumber'])
  wn_exp = df['wavenumber'].values

  n_water, k_water = get_water_nk(wn_exp, df_water_n_k)

    # check n of water
  plt.plot(wn_exp, n_water)
  plt.axhline(
    2.4 * np.sin(np.deg2rad(45)),
    linestyle="--"
  )
  plt.xlabel("Wavenumber (cm$^{-1}$)")
  plt.ylabel("Refractive index n")
  plt.show()

  qz = complex_qz(
    wn_exp,
    n_water,
    k_water
    )

  d_eff_cm = effective_depth_from_qz(qz)
  d_eff_um = d_eff_cm * 1e4


  # Check penetration depth
  plt.figure()
  plt.plot(wn_exp, d_eff_um)
  plt.xlabel("Wavenumber (cm$^{-1}$)")
  plt.ylabel("Effective depth (µm)")
  plt.title("Water-model ATR effective depth")
  plt.show()
  plt.figure()


  # save results
  export_array = np.stack((wn_exp, d_eff_cm)).T[0:]
  export_df = pd.DataFrame(export_array, columns=['wavenumber', 'effective penetration depth in cm'])

  try:
      with filedialog.asksaveasfile(mode='w', defaultextension=".csv") as file:
            export_df.to_csv(file.name, index=False)

  except TypeError:
        # if user cancels save, filedialog returns None rather than a file object, and the 'with' will raise an error
        print("The user cancelled save")

