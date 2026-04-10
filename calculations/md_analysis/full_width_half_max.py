import numpy as np
from scipy.signal import peak_widths, find_peaks
from scipy.ndimage import gaussian_filter1d

def calculate_fwhm(angles, distribution, sigma=2, prominence=0.01):
    """
    Calculate FWHM for all peaks in a bond angle distribution.
    
    Parameters
    ----------
    angles       : array of angle values (degrees)
    distribution : array of distribution values
    sigma        : gaussian smoothing width (default 2)
    prominence   : minimum peak prominence to detect (default 0.01)
    
    Returns
    -------
    results : list of dicts with peak position, FWHM, and crossing points
    """
    # smooth the distribution
    smoothed = gaussian_filter1d(distribution, sigma=sigma)

    # find peaks
    peaks, properties = find_peaks(smoothed, prominence=prominence)

    if len(peaks) == 0:
        print("No peaks found — try lowering prominence threshold")
        return []

    # calculate FWHM at each peak
    widths, width_heights, left_ips, right_ips = peak_widths(
        smoothed, peaks, rel_height=0.5
    )

    # convert from index units to degree units
    angle_spacing = angles[1] - angles[0]

    results = []
    for i, peak_idx in enumerate(peaks):
        result = {
            'peak_position' : angles[peak_idx],
            'peak_value'    : smoothed[peak_idx],
            'fwhm'          : widths[i] * angle_spacing,
            'left_crossing' : angles[0] + left_ips[i]  * angle_spacing,
            'right_crossing': angles[0] + right_ips[i] * angle_spacing,
            'half_max'      : width_heights[i],
        }
        results.append(result)

    return results


def print_fwhm_results(results, label=""):
    if label:
        print(f"\n{'='*40}")
        print(f"  {label}")
        print(f"{'='*40}")
    for i, r in enumerate(results):
        print(f"  Peak {i+1}:")
        print(f"    Position      : {r['peak_position']:.2f} degrees")
        print(f"    FWHM          : {r['fwhm']:.2f} degrees")
        print(f"    Left crossing : {r['left_crossing']:.2f} degrees")
        print(f"    Right crossing: {r['right_crossing']:.2f} degrees")


# ----------------------------------------------------------------
# example usage — replace with your actual data loading
# ----------------------------------------------------------------

# if loading from a file (two columns: angle, distribution)
# data = np.loadtxt("bond_angles.dat")
# angles = data[:, 0]
# hoh    = data[:, 1]
# oho    = data[:, 2]
# hoh2   = data[:, 3]

# placeholder: generate fake data for testing
angles = np.linspace(0, 180, 500)
hoh  = np.exp(-0.5 * ((angles - 104) / 5 )**2)   # H-O-H,  narrow ~104 deg
oho  = np.exp(-0.5 * ((angles - 165) / 15)**2)   # O-H--O, broad  ~165 deg
hoh2 = np.exp(-0.5 * ((angles - 110) / 20)**2)   # H--O--H, broad ~110 deg

filepath_h_o_h = '/Users/loganyamamoto/Desktop/Research/project_writeup/Vashishta Potential/figure_data/water_structure/ba_h-o-h_finetune_8973.dat'
data_h_o_h = np.loadtxt(filepath_h_o_h, delimiter=',', skiprows=1, usecols=range(9))
angles = data_h_o_h[:, 0]
hoh = data_h_o_h[:, 4]

filepath_o_h_o = '/Users/loganyamamoto/Desktop/Research/project_writeup/Vashishta Potential/figure_data/water_structure/ba_o-h--o_finetune_8973.dat'
data_o_h_o = np.loadtxt(filepath_o_h_o, delimiter=',', skiprows=1, usecols=range(9))
angles = data_o_h_o[:, 0]
oho = data_o_h_o[:, 5]

filepath_h_o_h2 = '/Users/loganyamamoto/Desktop/Research/project_writeup/Vashishta Potential/figure_data/water_structure/ba_h--o--h_finetune_8973.dat'
data_h_o_h2 = np.loadtxt(filepath_h_o_h2, delimiter=',', skiprows=1, usecols=range(9))
angles = data_h_o_h2[:, 0]
hoh2 = data_h_o_h2[:, 4]

# calculate and print FWHM for each bond angle type
for distribution, label in zip(
    [hoh, oho, hoh2],
    ["H-O-H", "O-H--O", "H--O--H"]
):
    results = calculate_fwhm(angles, distribution, sigma=2, prominence=0.01)
    print_fwhm_results(results, label=label)