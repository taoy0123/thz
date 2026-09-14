import pandas as pd

from tkinter.filedialog import askopenfilename
# from tkinter.messagebox import showinfo
from tkinter import filedialog


# Input and output file paths
input_file = askopenfilename()
scale_factor = float(input('Scaling factor:  '))

# Load the CSV file
df = pd.read_csv(input_file)

# Make a copy for modification
df_modified = df.copy()

# Select all columns starting from column B (index 1)
cols_to_modify = df.columns[1:]
# # Select all columns starting from column C (index 2)
# cols_to_modify = df.columns[2:]

# Multiply those columns by scaling factor (ignoring headers automatically)
df_modified.loc[:, cols_to_modify] = df_modified.loc[:, cols_to_modify] * scale_factor

# Export the result to a new CSV
try:
    with filedialog.asksaveasfile(mode='w', defaultextension=".csv") as file:
        df_modified.to_csv(file.name, index=False)

except TypeError:
    # if user cancels save, filedialog returns None rather than a file object, and the 'with' will raise an error
    print("The user cancelled save")