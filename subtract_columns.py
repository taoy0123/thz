import pandas as pd
from tkinter.filedialog import askopenfilename
from tkinter import filedialog


# File paths
file1 = 'C:/Users/Tao/Documents/lab_results/2025-05-14_Hank/water.csv'
file2 = askopenfilename()

# Read both CSVs
df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)

# Ensure both files have the same columns
if not df1.columns.equals(df2.columns):
    raise ValueError("The two CSV files must have identical column headers")

# Copy the first column unchanged
result = df1.iloc[:, :1].copy()

# Subtract from column B (index 1) onwards
for col in df1.columns[1:]:
    result[col] = df1[col] - df2[col]

# Save to CSV
try:
    with filedialog.asksaveasfile(mode='w', defaultextension=".csv") as file:
        result.to_csv(file.name, index=False)

except TypeError:
    # if user cancels save, filedialog returns None rather than a file object, and the 'with' will raise an error
    print("The user cancelled save")