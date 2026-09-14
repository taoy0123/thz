import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

from tkinter.filedialog import askopenfilename
from tkinter.messagebox import showinfo
from tkinter import filedialog

# Input files
showinfo(message='Select data file') 
file1 = askopenfilename()
showinfo(message='Select function file') 
file2 = askopenfilename()
# file2 = "da_Water_20250514.csv"


# Read CSVs
df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)

# Assumed format:
# col0 = sample name (ignored)
# col1 = frequency
# col2..end = spectra at different temperatures (headers are temperature labels)

# Verify matching frequency axis
if not np.allclose(df1.iloc[:,1].values, df2.iloc[:,1].values):
    raise ValueError("Frequency columns do not match between files")

freq = df1.iloc[:,1].values
temps = df1.columns[2:]   # temperature-dependent spectra headers

# Restrict fitting range: rows 23 to 594 (50 to 600 cm-1 inclusive)
row_start, row_end = 23, 594
fit_slice = slice(row_start, row_end + 1)

# Flatten all spectra pairs (spec1 vs spec2) into one long array for fitting
y1_all = []
y2_all = []

for temp in temps:
    y1_all.append(df1[temp].values[fit_slice])
    y2_all.append(df2[temp].values[fit_slice])

y1_all = np.concatenate(y1_all)
y2_all = np.concatenate(y2_all)

# # Determine the scalling factor by curve fitting
# # Define global scaling model
# def model(x, k):
#     return k * y2_all

# # Fit single scaling factor k_opt
# popt, _ = curve_fit(model, np.arange(len(y1_all)), y1_all)
# k_opt = popt[0]

# Determine the scalling factor by aligning the minima
k_opt = y1_all.min() / y2_all.min()

print(f"Global scaling factor k = {k_opt:.6f}")

# Compute residuals for each temperature
residuals = pd.DataFrame({"Frequency": freq})
for temp in temps:
    y1 = df1[temp].values
    y2 = df2[temp].values
    residuals[temp] = y1 - k_opt * y2

# Save residuals to CSV
try:
    with filedialog.asksaveasfile(mode='w', defaultextension=".csv") as file:
        residuals.to_csv(file.name, index=False)

except TypeError:
    # if user cancels save, filedialog returns None rather than a file object, and the 'with' will raise an error
    print("The user cancelled save")