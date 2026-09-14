import numpy as np
import pandas as pd
from tkinter.filedialog import askopenfilename
from tkinter import filedialog

# Load data
ref_df = pd.read_csv("da_Water_20250514.csv")
sample_df = pd.read_csv(askopenfilename())

wavenumber = ref_df.iloc[:,1].values
corrected_df = pd.DataFrame({"Wavenumber": wavenumber})

# Restrict to 50–600 cm-1
mask = (wavenumber >= 50) & (wavenumber <= 600)

for col in ref_df.columns[2:]:
    ref = ref_df[col].values
    sample = sample_df[col].values

    # Compute safe scaling factor: min(sample/ref) over masked region
    ratio = sample[mask] / ref[mask]
    ratio = ratio[ref[mask] > 0]  # avoid divide-by-zero
    coef = np.min(ratio)

    # Subtract baseline with restricted coefficient
    baseline = coef * ref
    corrected = sample - baseline

    corrected_df[col] = corrected

# Save result

try:
    with filedialog.asksaveasfile(mode='w', defaultextension=".csv") as file:
        corrected_df.to_csv(file.name, index=False)

except TypeError:
    # if user cancels save, filedialog returns None rather than a file object, and the 'with' will raise an error
    print("The user cancelled save")