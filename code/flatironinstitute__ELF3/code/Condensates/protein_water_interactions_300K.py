import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.ticker import MaxNLocator

def load_data(file_path):
    """ Load data from a single-column file without headers """
    return pd.read_csv(file_path, header=None)[0]

# Load data
data1 = load_data('water300K.dat')
data2 = load_data('protein300K.dat')

# Create a figure and a single subplot with a thinner aspect ratio
fig, ax1 = plt.subplots(figsize=(8, 3.2))  # Adjusted height to make it thinner

# Plot data1 on the primary y-axis
ax1.plot(data1, color='black')
ax1.set_ylabel('# Water Contacts', size=18)
ax1.tick_params(axis='y', labelsize=18)

# Create a secondary y-axis
ax2 = ax1.twinx()
ax2.plot(data2, color='red')
ax2.set_ylabel('# Protein Contacts', size=18)
ax2.tick_params(axis='y', labelsize=18)

# Set the x-axis labels
ax1.set_xlabel('Time ($\mu$s)', size=20)
ticks = np.linspace(0, len(data1)-1, num=7)  # Choose 7 ticks for example
labels = ['0.0', '0.5', '1', '1.5', '2', '2.5', '3.0']  # Labels for the ticks
ax1.set_xticks(ticks)
ax1.set_xticklabels(labels, fontsize=14)
ax1.set_xlim([0, len(data1)-1])

# Adjust tick label size
ax1.tick_params(axis='both', which='major', labelsize=18)
ax2.tick_params(axis='both', which='major', labelsize=18)

# Add text annotation
ax1.text(7500, 105, "7Q", size="26")

# Ensure the y-axis ticks are integers and reduce the number of tick labels
ax1.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
ax2.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))

# Adjust layout to prevent overlap and reduce space
plt.tight_layout()

# Show the plot
plt.show()
plt.savefig("manual_tick_labels.png", bbox_inches='tight')

