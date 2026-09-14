import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from tkinter.filedialog import askopenfilename
from tkinter.messagebox import showinfo
from tkinter import Tk, filedialog

# function file contains the spectrum to be scaled
showinfo(message='Select function file')  
sample_filename = askopenfilename()
# data file contains the spectrum to be subtracted from
showinfo(message='Select data file')
reference_filename = askopenfilename()

# Load data
ref = pd.read_csv(reference_filename)
target = pd.read_csv(sample_filename)

# Assume first column = wavenumber, second = ΔA_CO
wn_ref, delta_ref = ref.iloc[:, 0].values, ref.iloc[:, 1].values
wn_tar, delta_tar = target.iloc[:, 0].values, target.iloc[:, 1].values

# Interpolate both onto common grid
common_min = max(wn_ref.min(), wn_tar.min())
common_max = min(wn_ref.max(), wn_tar.max())
wn_common = np.linspace(common_min, common_max, min(len(wn_ref), len(wn_tar)))
interp_ref = np.interp(wn_common, wn_ref, delta_ref)
interp_tar = np.interp(wn_common, wn_tar, delta_tar)

# Select range for scaling
mask = (wn_common >= 50) & (wn_common <= 500)
wn_sel = wn_common[mask]
ref_sel = interp_ref[mask]
tar_sel = interp_tar[mask]

# Compute ratio where both are positive
valid = (ref_sel > 0) & (tar_sel > 0)
ratios = ref_sel[valid] / tar_sel[valid]

# IQR-based outlier removal
# -----------------------------
Q1 = np.quantile(ratios, 0.25)
Q3 = np.quantile(ratios, 0.75)
IQR = Q3 - Q1

# Keep only ratios inside the IQR window (1.5× rule)
lower_bound = Q1 - 1.5 * IQR
upper_bound = Q3 + 1.5 * IQR

ratio_mask = (ratios >= lower_bound) & (ratios <= upper_bound)
ratios_filtered = ratios[ratio_mask]

# If filtering removed everything, fall back to original ratios
if len(ratios_filtered) < 5:
    ratios_filtered = ratios

# Robust scale factor (now safe from low and high outliers)
scale_factor = np.quantile(ratios_filtered, 0.05)

print(f"Scaling factor of {sample_filename} is \n{scale_factor:.4f}")

# Scale full target
scaled_tar = delta_tar * scale_factor

# Compute difference on common grid (ref - scaled target)
scaled_interp_tar = np.interp(wn_common, wn_tar, scaled_tar)
diff = interp_ref - scaled_interp_tar

# Remove outliers (e.g., > 3σ)
mean_diff = np.mean(diff)
std_diff = np.std(diff)
mask_no_outliers = np.abs(diff - mean_diff) <= 3 * std_diff

wn_filtered = wn_common[mask_no_outliers]
diff_filtered = diff[mask_no_outliers]

# plot
range_min = 70
range_max = 500
plt.figure(figsize=(9, 5))
plt.subplot(2, 1, 1)
plt.plot(wn_tar, delta_tar, label="Target (original)", alpha=0.4)
plt.plot(wn_tar, scaled_tar, label=f"Target scaled ×{scale_factor:.4f}", lw=2)
plt.plot(wn_ref, delta_ref, label="Reference", lw=2, alpha=0.7)
plt.xlim(range_min, range_max)
plt.xlabel("Wavenumber (cm⁻¹)")
plt.ylabel("ΔA_CO")
plt.legend()
plt.title(f"Spectra ({range_min}-{range_max} cm⁻¹)")

plt.subplot(2, 1, 2)
plt.plot(wn_common, diff, label="Difference (raw)", color="gray", alpha=0.5)
plt.plot(wn_filtered, diff_filtered, label="Difference (filtered)", color="red", lw=1.5)
plt.xlim(range_min, range_max)
plt.xlabel("Wavenumber (cm⁻¹)")
plt.ylabel("ΔA_CO difference (Ref - Scaled)")
plt.legend()
plt.title("Difference spectrum (with and without outliers)")
plt.tight_layout()
plt.show()

# # Export results
# scaled_df = pd.DataFrame({"Wavenumber (cm^-1)": wn_tar, "ΔA_CO_scaled": scaled_tar})
# diff_df = pd.DataFrame({"Wavenumber (cm^-1)": wn_filtered, "Difference (Ref - Scaled)": diff_filtered})

# try:
#     with filedialog.asksaveasfile(mode='w', defaultextension=".csv") as file:
#         diff_df.to_csv(file.name, index=False)

# except TypeError:
#     # if user cancels save, filedialog returns None rather than a file object, and the 'with' will raise an error
#     print("The user cancelled save")