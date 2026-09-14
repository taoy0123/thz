'''The baseline is fit by scaling the reference with NNLS.

After subtraction, the corrected spectrum is processed so:

Small negative values (likely noise) are set to 0.

Real sharp dips (below threshold) remain visible as negative peaks.'''


import numpy as np
import pandas as pd
from scipy.optimize import nnls
from tkinter.filedialog import askopenfilename
from tkinter import filedialog

# ---- Load reference and sample ----
ref_file = "da_Water_20250514.csv"
sample_file = askopenfilename()

# Load as DataFrame
ref_df = pd.read_csv(ref_file)
sample_df = pd.read_csv(sample_file)

# Extract wavenumber
wavenumber = ref_df.iloc[:,1].values

# Prepare output DataFrame
corrected_df = pd.DataFrame()
corrected_df["Wavenumber"] = wavenumber

# Threshold for keeping sharp dips
threshold = -0.05

# Mask for the 50–600 cm-1 region
mask = (wavenumber >= 50) & (wavenumber <= 600)

# Iterate through spectra columns (from col 2 onward)
for col in ref_df.columns[2:]:
    ref = ref_df[col].values
    sample = sample_df[col].values

    # --- Fit baseline only in the 50–600 cm-1 region ---
    coef, _ = nnls(ref[mask, np.newaxis], sample[mask])
    baseline = coef * ref

    # --- Correct sample only in 50–600 region ---
    corrected = sample.copy()
    corrected[mask] = sample[mask] - baseline[mask]

    # # Apply positivity rule inside 50–600 cm-1
    # corrected[mask] = np.where(corrected[mask] < 0,
    #                            np.maximum(corrected[mask], threshold),
    #                            corrected[mask])

    # Add to output DataFrame
    corrected_df[col] = corrected

# ---- Save corrected spectra ----
try:
    with filedialog.asksaveasfile(mode='w', defaultextension=".csv") as file:
        corrected_df.to_csv(file.name, index=False)

except TypeError:
    # if user cancels save, filedialog returns None rather than a file object, and the 'with' will raise an error
    print("The user cancelled save")