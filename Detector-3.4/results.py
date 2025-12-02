import numpy as np

percentage_list = []

# === percentage list and average ===
def add_to_list(percentage, out_handle, final=False):
    if not final:
        percentage_list.append(round(percentage, 3))
        if percentage >= 95.0:
            out_handle.write(f"This sample is {percentage:.2f}% human, indicating it is a human sample.\n\n")
        else:
            out_handle.write(f"This sample is {percentage:.2f}% human, indicating it is not a human sample.\n\n")
    else:
        if percentage_list:
            percent_average = np.mean(percentage_list)
            out_handle.write(f"\nOn average, this sample is {percent_average:.2f}% human!\n")
        else:
            out_handle.write("\nNo valid alignments were performed.\n")
