"""
Visualization of Flory-Huggins Prewetting Transition
Heatmap of largest_droplet_size averaged over final 20% of frames
for varying temperature T and concentration α
"""

import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import glob
import sys

# Configuration
DATA_DIR = Path("/home/tohme/.julia/dev/Droplet/data")
TEMPERATURES = [0.56, 0.58, 0.60, 0.62, 0.64, 0.66]
NCPU = 32  # Number of CPU runs per temperature (32 concentrations)

def load_and_aggregate_data(date_folders, data_key="largest_droplet_size"):
    """
    Load JSON data from all date/temperature combinations and aggregate.
    
    Parameters:
        date_folders: list of date folder names to aggregate
    
    Returns:
        heatmap_data: (n_temps, n_cpus) array of averaged largest_droplet_size
        concentrations: list of concentration values from filenames
    """
    
    # Dictionary to store aggregated data: {temp: {cpu_idx: [values]}}
    aggregated = {t: {i: [] for i in range(1, NCPU + 1)} for t in TEMPERATURES}
    
    print(f"Scanning data directories: {date_folders}")
    print(f"Temperatures: {TEMPERATURES}")
    print(f"Expected CPU runs per temperature: {NCPU}")
    
    count_by_temp = {t: 0 for t in TEMPERATURES}
    
    # Iterate through date folders
    for date_folder in date_folders:
        date_path = DATA_DIR / date_folder
        if not date_path.exists():
            print(f"Warning: {date_path} does not exist")
            continue
            
        # Iterate through temperature folders
        for temp in TEMPERATURES:
            temp_folder_pattern = f"flory_huggins_prewetting_42_T_{temp:.2f}"
            temp_path = date_path / temp_folder_pattern
            
            if not temp_path.exists():
                print(f"Warning: {temp_path} does not exist")
                continue
            
            # Load data from each CPU output file
            for cpu_idx in range(1, NCPU + 1):
                json_file = temp_path / f"output_cpu{cpu_idx}.json"
                
                if not json_file.exists():
                    print(f"Warning: {json_file} does not exist")
                    continue
                
                try:
                    with open(json_file, 'r') as f:
                        data = json.load(f)
                    
                    # Data is nested with two numeric keys at top level
                    first_key = list(data.keys())[0]
                    second_key = list(data[first_key].keys())[0]
                    sim_data = data[first_key][second_key]
                    
                    # Extract necessary fields
                    largest_droplet_size = sim_data[data_key]
                    nframes = sim_data['nframes']
                    concentration = sim_data['concentration']
                    
                    # Calculate threshold for final 20% of frames
                    threshold_frame = int(0.8 * nframes)
                    
                    # Average over final 20% of frames
                    avg_size = np.mean(largest_droplet_size[threshold_frame:])
                    
                    # Store aggregated value
                    aggregated[temp][cpu_idx].append({
                        'avg_size': avg_size,
                        'concentration': concentration,
                        'date': date_folder
                    })
                    
                    count_by_temp[temp] += 1
                    
                except Exception as e:
                    print(f"Error processing {json_file}: {e}")
                    continue
    
    print(f"\nData loading summary:")
    for temp, count in count_by_temp.items():
        print(f"  T={temp}: {count} files loaded")
    
    # Now organize data into heatmap format
    # Shape: (n_temps, n_cpus) where each element is the average over dates
    heatmap_data = np.zeros((len(TEMPERATURES), NCPU))
    concentrations_by_cpu = {}
    
    for temp_idx, temp in enumerate(TEMPERATURES):
        for cpu_idx in range(1, NCPU + 1):
            values = aggregated[temp][cpu_idx]
            if values:
                # Average across dates
                avg_across_dates = np.mean([v['avg_size'] for v in values])
                heatmap_data[temp_idx, cpu_idx - 1] = avg_across_dates
                
                # Store concentration (should be same for all dates)
                if cpu_idx not in concentrations_by_cpu:
                    concentrations_by_cpu[cpu_idx] = values[0]['concentration']
            else:
                heatmap_data[temp_idx, cpu_idx - 1] = np.nan
    
    # Extract concentrations in order
    concentrations = [concentrations_by_cpu[i] for i in range(1, NCPU + 1)]
    
    # Sort by concentration to ensure proper ordering
    sorted_indices = np.argsort(concentrations)
    heatmap_data_sorted = heatmap_data[:, sorted_indices]
    concentrations_sorted = np.array(concentrations)[sorted_indices]
    
    return heatmap_data_sorted, concentrations_sorted

def compute_staircase_boundary(heatmap_data, concentrations, temperatures, threshold):
    """
    Compute the staircase boundary where values transition above/below threshold.
    Returns points that form a true staircase aligned with grid cell boundaries.
    
    Parameters:
        heatmap_data: (n_temps, n_concs) array
        concentrations: array of concentration values
        temperatures: list of temperature values
        threshold: value to threshold at
    
    Returns:
        staircase_concs: concentration values for staircase points
        staircase_temps: temperature values for staircase points
    """
    
    # For each temperature, find the first concentration where value exceeds threshold
    crossing_concs = {}  # temp_idx -> conc_idx
    
    for temp_idx, temp in enumerate(temperatures):
        row = heatmap_data[temp_idx, :]
        above_threshold = row > threshold
        
        if np.any(above_threshold):
            crossing_concs[temp_idx] = np.argmax(above_threshold)
    
    if not crossing_concs:
        return np.array([]), np.array([])
    
    # Build staircase by connecting crossings at grid boundaries
    final_concs = []
    final_temps = []
    
    sorted_temp_indices = sorted(crossing_concs.keys())
    
    for i, temp_idx in enumerate(sorted_temp_indices):
        conc_idx = crossing_concs[temp_idx]
        temp = temperatures[temp_idx]
        
        # Cell boundary between conc_idx-1 and conc_idx
        if conc_idx > 0:
            conc_boundary = (concentrations[conc_idx - 1] + concentrations[conc_idx]) / 2
        else:
            conc_boundary = concentrations[0]
        
        if i == 0:
            # First point
            final_concs.append(conc_boundary)
            final_temps.append(temp)
        else:
            # Previous temperature index
            prev_temp_idx = sorted_temp_indices[i - 1]
            prev_temp = temperatures[prev_temp_idx]
            prev_conc_idx = crossing_concs[prev_temp_idx]
            
            if prev_conc_idx > 0:
                prev_conc_boundary = (concentrations[prev_conc_idx - 1] + concentrations[prev_conc_idx]) / 2
            else:
                prev_conc_boundary = concentrations[0]
            
            # Temperature boundary between prev_temp_idx and temp_idx
            temp_boundary = (prev_temp + temp) / 2
            
            # Horizontal segment at temp_boundary
            final_concs.append(conc_boundary)
            final_temps.append(temp_boundary)
            
            # Vertical segment
            final_concs.append(conc_boundary)
            final_temps.append(temp)
    
    return np.array(final_concs), np.array(final_temps)

def plot_heatmap(heatmap_data, concentrations, temperatures, date_str, file_name, 
                staircase_concs=None, staircase_temps=None):
    """
    Create and save heatmap visualization.
    
    Parameters:
        heatmap_data: (n_temps, n_concs) array
        concentrations: array of concentration values
        temperatures: list of temperature values
        date_str: date folder name to append to filename
        file_name: output filename
        staircase_concs: optional concentration values for staircase overlay
        staircase_temps: optional temperature values for staircase overlay
    """
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
    })
    
    fig, ax = plt.subplots(figsize=(14, 6))
    
    # Create heatmap
    im = ax.imshow(heatmap_data, aspect='auto', cmap='viridis', origin='lower',
                   extent=[concentrations[0], concentrations[-1], 
                          temperatures[0], temperatures[-1]],
                   interpolation='nearest')
    
    # Overlay staircase if provided
    if staircase_concs is not None and staircase_temps is not None:
        ax.plot(staircase_concs, staircase_temps, 'k-', linewidth=2.5, label='Control Threshold')
        ax.legend(loc='upper left', fontsize=10)
    
    # Labels and title
    ax.set_xlabel('Concentration (α)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Temperature (T)', fontsize=12, fontweight='bold')
    ax.set_title('Flory-Huggins Prewetting Transition\nAverage Largest Droplet Size (Final 20% of Frames)',
                fontsize=13, fontweight='bold', pad=20)
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax, label='Largest Droplet Size (lattice units)')
    
    # Set tick positions
    ax.set_yticks(temperatures)
    ax.set_xticks(np.arange(0.005, 0.035, step=0.005))
    
    # Formatting
    ax.grid(False)
    fig.tight_layout()
    
    # Save figure
    output_path = Path(f"/home/tohme/.julia/dev/Droplet/media/figures/{date_str}/{file_name}.pdf")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nHeatmap saved to: {output_path}")
    
    plt.close()

def main():
    print("=" * 70)
    print("Flory-Huggins Prewetting Transition Visualization")
    print("=" * 70)
    
    # Get arguments from command-line
    if len(sys.argv) < 4:
        print("Usage: python prewetting_heatmap.py <date_folder> <file_name> <data_key> [control_folder] [threshold]")
        print("Example: python prewetting_heatmap.py 20260131 output largest_droplet_size")
        print("Example with control: python prewetting_heatmap.py 20260131 output largest_droplet_size 20260130 3.5")
        sys.exit(1)
    
    date_folders = [sys.argv[1]]
    file_name = sys.argv[2]
    data_key = sys.argv[3]
    control_folder = sys.argv[4] if len(sys.argv) > 4 else None
    threshold = float(sys.argv[5]) if len(sys.argv) > 5 else None
    
    # Load and aggregate main data
    heatmap_data, concentrations = load_and_aggregate_data(date_folders, data_key)
    
    # Load control data if provided
    staircase_concs = None
    staircase_temps = None
    
    if control_folder is not None and threshold is not None:
        print(f"\nLoading control data from: {control_folder}")
        print(f"Threshold value: {threshold}")
        control_data, control_concs = load_and_aggregate_data([control_folder], data_key)
        staircase_concs, staircase_temps = compute_staircase_boundary(
            control_data, control_concs, TEMPERATURES, threshold
        )
        print(f"Staircase boundary computed with {len(staircase_temps)} points")
    
    # Print statistics
    print(f"\nData statistics:")
    print(f"  Heatmap shape: {heatmap_data.shape}")
    print(f"  Temperature range: {TEMPERATURES[0]:.2f} - {TEMPERATURES[-1]:.2f}")
    print(f"  Concentration range: {concentrations[0]:.6f} - {concentrations[-1]:.6f}")
    print(f"  Mean droplet size: {np.nanmean(heatmap_data):.2f}")
    print(f"  Std droplet size: {np.nanstd(heatmap_data):.2f}")
    print(f"  Min droplet size: {np.nanmin(heatmap_data):.2f}")
    print(f"  Max droplet size: {np.nanmax(heatmap_data):.2f}")
    
    # Create and save heatmap
    plot_heatmap(heatmap_data, concentrations, TEMPERATURES, date_folders[0], file_name,
                staircase_concs, staircase_temps)
    
    print("\n" + "=" * 70)
    print("Visualization complete!")
    print("=" * 70)

if __name__ == "__main__":
    main()
