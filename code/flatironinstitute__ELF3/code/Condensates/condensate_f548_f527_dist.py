import pandas as pd
import matplotlib.pyplot as plt

# Function to read data from a file
def read_data(filename):
    # Assuming the data is in a single column without header
    return pd.read_csv(filename, header=None).iloc[:, 0]

# File paths (replace with your actual file paths)
file1 = 'dist_all.dat'
file2 = 'dist_all340K.dat'

# Read the data
data1 = read_data(file1)
data2 = read_data(file2)

# Plotting
plt.figure(figsize=(10, 6))

# Create histograms
plt.hist(data1, bins=150, alpha=0.5, label='290K', color="black")
plt.hist(data2, bins=150, alpha=0.5, label='340K', color="red")

# Setting the labels, title, and legend
plt.xlabel('F458-F527 Distance (nm)', fontsize=20)
plt.ylabel('# Conformations', fontsize=20)
#plt.title('Histogram of Distances', fontsize=20)
plt.legend(fontsize=20)

# Adjusting tick size
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

# Saving the figure (optional)
plt.savefig('histograms2.png', dpi=300,bbox_inches="tight")

# Show plot
plt.show()
