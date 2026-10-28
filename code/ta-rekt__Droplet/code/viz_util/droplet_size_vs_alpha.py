"""
Plot maximum droplet size vs nucleation site size (alpha)
Averages largest_droplet_size over the final 20% of frames for each alpha
"""

import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# Configuration
DATA_DIR = Path("/home/tohme/.julia/dev/Droplet/data")

def load_and_process_data(date_folder, temp_folder):
    """
    Load JSON data from prewetting_size_scaling simulations and process.
    
    Parameters:
        date_folder: date folder name (e.g., "20260131")
        temperature: temperature value (e.g., 0.66)
    
    Returns:
        alphas: array of nucleation site sizes
        avg_max_sizes: array of average largest droplet sizes (final 20% of frames)
    """
    
    date_path = DATA_DIR / date_folder
    temp_path = date_path / temp_folder
    
    if not temp_path.exists():
        print(f"Error: {temp_path} does not exist")
        return None, None
    
    alphas = []
    avg_max_sizes = []
    
    # Load all output files
    json_files = sorted(temp_path.glob("output_cpu*.json"))
    
    if not json_files:
        print(f"Error: No JSON files found in {temp_path}")
        return None, None
    
    print(f"Found {len(json_files)} output files")
    
    for json_file in json_files:
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
            
            # Data is nested with two numeric keys at top level
            first_key = list(data.keys())[0]
            second_key = list(data[first_key].keys())[0]
            sim_data = data[first_key][second_key]
            
            # Extract necessary fields
            largest_droplet_size = sim_data['largest_droplet_size']
            nframes = sim_data['nframes']
            alpha = sim_data['site_length']
            
            # Calculate threshold for final 20% of frames
            threshold_frame = int(0.8 * nframes)
            
            # Average over final 20% of frames
            avg_size = np.mean(largest_droplet_size[threshold_frame:])
            
            alphas.append(alpha)
            avg_max_sizes.append(avg_size)
            
            print(f"  Loaded {json_file.name}: alpha={alpha}, avg_size={avg_size:.2f}")
            
        except Exception as e:
            print(f"Error processing {json_file}: {e}")
            continue
    
    # Sort by alpha
    sorted_indices = np.argsort(alphas)
    alphas = np.array(alphas)[sorted_indices]
    avg_max_sizes = np.array(avg_max_sizes)[sorted_indices]
    
    return alphas, avg_max_sizes

def plot_droplet_size_vs_alpha(alphas, avg_max_sizes, temperature, conc, date_folder):
    """
    Create and save plot of droplet size vs nucleation site size.
    
    Parameters:
        alphas: array of nucleation site sizes
        avg_max_sizes: array of average largest droplet sizes
        temperature: temperature value
        date_folder: date folder name for output filename
    """
    
    fig, ax = plt.subplots(figsize=(7, 5))
    
    # Plot
    ax.plot(alphas, avg_max_sizes, 'o-', linewidth=2, markersize=8, color='steelblue')
    ax.set_ylim(0, 2*np.mean(avg_max_sizes))
    
    # Labels and title
    ax.set_xlabel('Nucleation Site Size (α)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Average Largest Droplet Size (Final 20% of Frames)', fontsize=12, fontweight='bold')
    ax.set_title(f'Droplet Size vs Nucleation Site Size\nT = {temperature:.2f}\nc = {conc:.3f}',
                fontsize=13, fontweight='bold', pad=20)
    
    # Grid
    ax.grid(True, alpha=0.3)
    
    # Formatting
    fig.tight_layout()
    
    # Save figure
    output_dir = Path(f"/home/tohme/.julia/dev/Droplet/media/figures/{date_folder}")
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"droplet_size_vs_alpha_c_{conc:.3f}_T_{temperature:.2f}.png"
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nPlot saved to: {output_path}")
    
    plt.close()

def main():
    print("=" * 70)
    print("Droplet Size vs Nucleation Site Size")
    print("=" * 70)
    
    # Get date folder and temperature from command-line arguments
    if len(sys.argv) < 3:
        print("Usage: python droplet_size_vs_alpha.py <date_folder> <temperature>")
        print("Example: python droplet_size_vs_alpha.py 20260131 0.66")
        sys.exit(1)
    
    date_folder = sys.argv[1]
    temp_folder = sys.argv[2]
    temperature = float(temp_folder[-4:])
    conc = float(temp_folder[-12:-7])
    
    # Load and process data
    alphas, avg_max_sizes = load_and_process_data(date_folder, temp_folder)
    
    if alphas is None or avg_max_sizes is None:
        print("Failed to load data")
        sys.exit(1)
    
    # Print statistics
    print(f"\nData statistics:")
    print(f"  Number of alpha values: {len(alphas)}")
    print(f"  Alpha range: {alphas[0]:.6f} - {alphas[-1]:.6f}")
    print(f"  Average droplet size range: {np.min(avg_max_sizes):.2f} - {np.max(avg_max_sizes):.2f}")
    print(f"  Mean average droplet size: {np.mean(avg_max_sizes):.2f}")
    
    # Create and save plot
    plot_droplet_size_vs_alpha(alphas, avg_max_sizes, temperature, conc, date_folder)
    
    print("\n" + "=" * 70)
    print("Plot complete!")
    print("=" * 70)

if __name__ == "__main__":
    main()
