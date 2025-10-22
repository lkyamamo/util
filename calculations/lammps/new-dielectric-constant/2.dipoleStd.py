import sys
import csv
from collections import deque

file1 = sys.argv[1]
T = float(sys.argv[2])
la = float(sys.argv[3])
lb = float(sys.argv[4])
lc = float(sys.argv[5])

# Window size for averaging
WINDOW_SIZE = 10000

# Bin size for binned averaging
BIN_SIZE = 100

# Pre-compute constants for better performance
WINDOW_SIZE_FLOAT = float(WINDOW_SIZE)
DEBYE_CONVERSION_FACTOR = 4.8032**2  # Pre-compute the conversion factor

# Choose averaging method: 'windowed', 'cumulative', 'hybrid', or 'binned'
# 'windowed': Only calculate averages when window is full (10,000 frames)
# 'cumulative': Calculate averages over all frames up to current point
# 'hybrid': Use cumulative averaging until window is full, then switch to windowed
# 'binned': Use standalone binned averaging (every 100 frames)
AVERAGING_METHOD = 'binned'  # Change to 'cumulative', 'hybrid', or 'binned' if desired

# Timing configuration for large datasets
TIMING_CHECKPOINT_INTERVAL = 50000  # Progress updates every N frames (adjust for dataset size)

def calculate_windowed_averages(dipole_window, running_sum_msq=0.0, running_sum_mx=0.0, running_sum_my=0.0, running_sum_mz=0.0):
    """
    Calculate averages over the sliding window using running sums.
    Returns None if window is not full yet.
    """
    if len(dipole_window) < WINDOW_SIZE:
        return None
    
    # Use running sums directly (much more efficient)
    window_msq = running_sum_msq / (3*WINDOW_SIZE_FLOAT)
    window_mx = running_sum_mx / WINDOW_SIZE_FLOAT
    window_my = running_sum_my / WINDOW_SIZE_FLOAT
    window_mz = running_sum_mz / WINDOW_SIZE_FLOAT
    
    window_mavgsq = (window_mx*window_mx + window_my*window_my + window_mz*window_mz)/3.0
    window_var = (window_msq - window_mavgsq)  # Combine operations
    
    return {
        'window_msq': window_msq,
        'window_mx': window_mx,
        'window_my': window_my,
        'window_mz': window_mz,
        'window_mavgsq': window_mavgsq,
        'window_var': window_var
    }

def calculate_cumulative_averages(dipole_window, running_sum_msq=0.0, running_sum_mx=0.0, running_sum_my=0.0, running_sum_mz=0.0):
    """
    Calculate cumulative averages over all frames up to current point using running sums.
    Returns averages even if window is not full.
    """
    if len(dipole_window) == 0:
        return None
    
    # Use running sums directly (much more efficient)
    frame_count = float(len(dipole_window))  # Convert once
    cumul_msq = running_sum_msq / (3*frame_count)
    cumul_mx = running_sum_mx / frame_count
    cumul_my = running_sum_my / frame_count
    cumul_mz = running_sum_mz / frame_count
    
    cumul_mavgsq = (cumul_mx*cumul_mx + cumul_my*cumul_my + cumul_mz*cumul_mz)/3.0
    cumul_var = (cumul_msq - cumul_mavgsq)  # Combine operations
    
    return {
        'cumul_msq': cumul_msq,
        'cumul_mx': cumul_mx,
        'cumul_my': cumul_my,
        'cumul_mz': cumul_mz,
        'cumul_mavgsq': cumul_mavgsq,
        'cumul_var': cumul_var
    }

def calculate_hybrid_averages(dipole_window, running_sum_msq=0.0, running_sum_mx=0.0, running_sum_my=0.0, running_sum_mz=0.0):
    """
    Hybrid approach: Use cumulative averaging until window is full,
    then switch to windowed averaging using running sums.
    """
    if len(dipole_window) == 0:
        return None
    
    # Once deque reaches WINDOW_SIZE, it never exceeds it due to maxlen
    # So we only need to check if it's exactly WINDOW_SIZE
    if len(dipole_window) == WINDOW_SIZE:
        # Use windowed averaging (deque with maxlen automatically maintains WINDOW_SIZE)
        hybrid_msq = running_sum_msq / (3*WINDOW_SIZE_FLOAT)
        hybrid_mx = running_sum_mx / WINDOW_SIZE_FLOAT
        hybrid_my = running_sum_my / WINDOW_SIZE_FLOAT
        hybrid_mz = running_sum_mz / WINDOW_SIZE_FLOAT
    else:
        # Use cumulative averaging until window is full
        frame_count = float(len(dipole_window))  # Convert once
        hybrid_msq = running_sum_msq / (3*frame_count)
        hybrid_mx = running_sum_mx / frame_count
        hybrid_my = running_sum_my / frame_count
        hybrid_mz = running_sum_mz / frame_count
    
    hybrid_mavgsq = (hybrid_mx*hybrid_mx + hybrid_my*hybrid_my + hybrid_mz*hybrid_mz)/3.0
    hybrid_var = (hybrid_msq - hybrid_mavgsq)  # Combine operations
    
    return {
        'hybrid_msq': hybrid_msq,
        'hybrid_mx': hybrid_mx,
        'hybrid_my': hybrid_my,
        'hybrid_mz': hybrid_mz,
        'hybrid_mavgsq': hybrid_mavgsq,
        'hybrid_var': hybrid_var
    }


def process_dipole_data(file1, averaging_function, data_key_prefix):
    """
    Process dipole moment data using the specified averaging function.
    
    Args:
        file1: Input file path
        averaging_function: Function to calculate averages (windowed, cumulative, or hybrid)
        data_key_prefix: Prefix for data keys in the returned dictionary
    
    Returns:
        tuple: (timestep_data, mavg, msq, counter, minitial, mxx, myy, mzz)
    """
    # Initialize sliding window
    dipole_window = deque(maxlen=WINDOW_SIZE)
    
    # Running sums for efficient window calculations
    running_sum_msq = 0.0
    running_sum_mx = 0.0
    running_sum_my = 0.0
    running_sum_mz = 0.0
    
    # Data collection for CSV output
    timestep_data = []
    
    # Variables for final calculations (over all data)
    mavg = [0.0, 0.0, 0.0]
    msq = 0.0
    counter = 0.0
    minitial = [0.0, 0.0, 0.0]
    mxx, myy, mzz = 0, 0, 0
    
    with open(file1, "r") as f1:
        for line in f1:
            line = line.strip().split()
            
            step, m1, m2, m3 = int(line[0]), float(line[1]), float(line[2]), float(line[3])
            
            # Track initial values
            if counter == 0:
                minitial[0], minitial[1], minitial[2] = m1, m2, m3
            
            # Add current dipole to sliding window and update running sums
            old_element = None
            if len(dipole_window) == WINDOW_SIZE:
                # Window is full, so the oldest element will be removed
                old_element = dipole_window[0]  # Get the element that will be removed
            
            dipole_window.append((m1, m2, m3))
            
            # Update running sums
            # Add new element
            running_sum_msq += m1*m1 + m2*m2 + m3*m3
            running_sum_mx += m1
            running_sum_my += m2
            running_sum_mz += m3
            
            # Remove old element if window was full
            if old_element is not None:
                old_m1, old_m2, old_m3 = old_element
                running_sum_msq -= old_m1*old_m1 + old_m2*old_m2 + old_m3*old_m3
                running_sum_mx -= old_m1
                running_sum_my -= old_m2
                running_sum_mz -= old_m3
            
            # Calculate averages using the provided function with running sums
            avg_data = averaging_function(dipole_window, running_sum_msq, running_sum_mx, running_sum_my, running_sum_mz)
            
            # Store data for CSV output
            if avg_data is not None:
                # Pre-compute Debye conversions to avoid repeated calculations
                mavgsq_debye = avg_data[f'{data_key_prefix}_mavgsq'] * DEBYE_CONVERSION_FACTOR
                msq_debye = avg_data[f'{data_key_prefix}_msq'] * DEBYE_CONVERSION_FACTOR
                var_debye = avg_data[f'{data_key_prefix}_var'] * DEBYE_CONVERSION_FACTOR
                
                timestep_data.append({
                    'timestep': step,
                    'mavgsq_t((eA)^2)': avg_data[f'{data_key_prefix}_mavgsq'],
                    'mavgsq_t(D^2)': mavgsq_debye,
                    'msqavg_t((eA)^2)': avg_data[f'{data_key_prefix}_msq'],
                    'msqavg_t(D^2)': msq_debye,
                    'var_t((eA)^2)': avg_data[f'{data_key_prefix}_var'],
                    'var_t(D^2)': var_debye
                })
            
            # Continue accumulating for final calculations (over all data)
            # Pre-compute squares to avoid repeated multiplication
            m1_sq = m1 * m1
            m2_sq = m2 * m2
            m3_sq = m3 * m3
            
            msq += m1_sq + m2_sq + m3_sq
            mxx += m1_sq
            myy += m2_sq
            mzz += m3_sq
            
            mavg[0] += m1
            mavg[1] += m2
            mavg[2] += m3
            
            counter += 1.0
            
            # Timing checkpoint: Simple progress indicator
            if counter % TIMING_CHECKPOINT_INTERVAL == 0:
                print(f"Processing frame {int(counter):,}")
    
    # Timing checkpoint: Processing complete
    print(f"Processing complete: {int(counter):,} frames")
    
    return timestep_data, mavg, msq, counter, minitial, mxx, myy, mzz

def process_dipole_data_binned(file1, bin_size=BIN_SIZE):
    """
    Process dipole moment data using binned averaging with running sums.
    Creates bins of specified size and averages all values within each bin.
    Uses efficient running sums instead of storing all bin data in memory.
    
    Args:
        file1: Input file path
        bin_size: Number of frames per bin (default: BIN_SIZE)
    
    Returns:
        tuple: (binned_data, mavg, msq, counter, mxx, myy, mzz)
    """
    # Data collection for binned output
    binned_data = []
    
    # Variables for final calculations (over all data)
    mavg = [0.0, 0.0, 0.0]
    msq = 0.0
    counter = 0.0
    mxx, myy, mzz = 0, 0, 0
    
    # Running sums for current bin (much more efficient)
    bin_sum_msq = 0.0
    bin_sum_mx = 0.0
    bin_sum_my = 0.0
    bin_sum_mz = 0.0
    bin_counter = 0
    last_step = 0  # Track the last timestep for the current bin
    
    with open(file1, "r") as f1:
        for line in f1:
            line = line.strip().split()
            
            step, m1, m2, m3 = int(line[0]), float(line[1]), float(line[2]), float(line[3])
            
            # Add to current bin running sums
            bin_sum_msq += m1*m1 + m2*m2 + m3*m3
            bin_sum_mx += m1
            bin_sum_my += m2
            bin_sum_mz += m3
            bin_counter += 1
            last_step = step  # Keep track of the latest timestep in current bin
            
            # When bin is full, calculate averages
            if bin_counter == bin_size:
                
                # Calculate averages for this bin
                bin_size_float = float(bin_size)
                binned_msq = bin_sum_msq / (3*bin_size_float)
                binned_mx = bin_sum_mx / bin_size_float
                binned_my = bin_sum_my / bin_size_float
                binned_mz = bin_sum_mz / bin_size_float
                
                binned_mavgsq = (binned_mx*binned_mx + binned_my*binned_my + binned_mz*binned_mz)/3.0
                binned_var = (binned_msq - binned_mavgsq) 
                
                # Use the last timestep in the bin as the representative timestep
                representative_step = last_step
                
                # Pre-compute Debye conversions
                mavgsq_debye = binned_mavgsq * DEBYE_CONVERSION_FACTOR
                msq_debye = binned_msq * DEBYE_CONVERSION_FACTOR
                var_debye = binned_var * DEBYE_CONVERSION_FACTOR
                
                binned_data.append({
                    'timestep': representative_step,
                    'mavgsq_t((eA)^2)': binned_mavgsq,
                    'mavgsq_t(D^2)': mavgsq_debye,
                    'msqavg_t((eA)^2)': binned_msq,
                    'msqavg_t(D^2)': msq_debye,
                    'var_t((eA)^2)': binned_var,
                    'var_t(D^2)': var_debye
                })
                
                # Reset bin running sums for next iteration
                bin_sum_msq = 0.0
                bin_sum_mx = 0.0
                bin_sum_my = 0.0
                bin_sum_mz = 0.0
                bin_counter = 0
            
            # Continue accumulating for final calculations (over all data)
            # Pre-compute squares to avoid repeated multiplication
            m1_sq = m1 * m1
            m2_sq = m2 * m2
            m3_sq = m3 * m3
            
            msq += m1_sq + m2_sq + m3_sq
            mxx += m1_sq
            myy += m2_sq
            mzz += m3_sq
            
            mavg[0] += m1
            mavg[1] += m2
            mavg[2] += m3
            
            counter += 1.0
            
            # Timing checkpoint: Simple progress indicator
            if counter % TIMING_CHECKPOINT_INTERVAL == 0:
                print(f"Processing frame {int(counter):,}")
    
    # Timing checkpoint: Processing complete
    print(f"Processing complete: {int(counter):,} frames")
    
    return binned_data, mavg, msq, counter, mxx, myy, mzz

# Select the appropriate averaging function based on the method
# Check binned method first
if AVERAGING_METHOD == 'binned':
    # Binned method uses standalone function, no averaging_function needed
    data_key_prefix = 'binned'
    timestep_data, mavg, msq, counter, mxx, myy, mzz = process_dipole_data_binned(file1)
else:
    if AVERAGING_METHOD == 'windowed':
        averaging_function = calculate_windowed_averages
        data_key_prefix = 'window'
    elif AVERAGING_METHOD == 'cumulative':
        averaging_function = calculate_cumulative_averages
        data_key_prefix = 'cumul'
    elif AVERAGING_METHOD == 'hybrid':
        averaging_function = calculate_hybrid_averages
        data_key_prefix = 'hybrid'
    else:
        raise ValueError(f"Unknown averaging method: {AVERAGING_METHOD}")

    timestep_data, mavg, msq, counter, minitial, mxx, myy, mzz = process_dipole_data(
        file1, averaging_function, data_key_prefix
    )



msq /= float(counter)

mxx /= float(counter)
myy /= float(counter)
mzz /= float(counter)

mavg[0] /= float(counter)
mavg[1] /= float(counter)
mavg[2] /= float(counter)

mdev = msq - (mavg[0]*mavg[0] + mavg[1]*mavg[1] + mavg[2]*mavg[2])
mx_dev = mxx - mavg[0]**2
my_dev = myy - mavg[1]**2
mz_dev = mzz - mavg[2]**2
mdev = mdev / 3

print("%14.6f %14.6f %14.6f" %(mavg[0], mavg[1],mavg[2]))
"""
with open(file1,"r") as f1:
	for line in f1.readlines():
		line = line.strip().split()
		m1,m2,m3 = float(line[0]), float(line[1]), float(line[2])
		m1,m2,m3 = m1-mavg[0],m2-mavg[1],m3-mavg[2]
		msq += (m1*m1) + (m2*m2) + (m3*m3)
"""

# units of K A^3/(e^2 m^2)
prefactor = (2.56*(10**7))/(1.38*8.85*T*la*lb*lc)

# Write timestep data to CSV
output_filename = file1.replace('.txt', f'_timestep_data_{AVERAGING_METHOD}.csv')
with open(output_filename, 'w', newline='') as csvfile:
    fieldnames = ['timestep', 'mavgsq_t((eA)^2)', 'mavgsq_t(D^2)', 'msqavg_t((eA)^2)', 'msqavg_t(D^2)', 'var_t((eA)^2)', 'var_t(D^2)']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    
    writer.writeheader()
    for row in timestep_data:
        writer.writerow(row)

print(f"Timestep data saved to: {output_filename}")

print("x_dev = %12.6f, y_dev = %12.6f, z_dev = %12.6f,deviation = %14.6f" %(mx_dev, my_dev, mz_dev, mdev))
print("eps_x = %12.6f, eps_y = %12.6f, eps_z = %12.6f, eps_total = %12.6f" %(prefactor * mx_dev, prefactor * my_dev, prefactor * mz_dev, prefactor * mdev))
