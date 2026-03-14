import numpy as np
import os
import h5py
import matplotlib.pyplot as plt
from scipy.ndimage import convolve1d
from scipy.signal import welch
import sys
sys.path.append("back_to_master1220/Moist_Dycore/IdealizeSpetral.jl/exp/HSt42/")
from EOF import EOF

# ==========================================
# 1. Configuration & Experiment Setup
# ==========================================
Z_DIM = 20  
Y_DIM = 32  
X_DIM = 128 
DAYS_PER_FILE = 25  
STEPS_PER_FILE = 100
EARTH_RADIUS = 6371000.0

start_days = range(500, 20000, 25)
total_days = len(start_days) * DAYS_PER_FILE

print(f"Total files per experiment: {len(start_days)}")
print(f"Total daily steps: {total_days}")

# Define all experiments
experiments = []
PR_values = [10, 20, 30, 40, 50]
pr_colors = np.array([[112, 115, 115], [182, 203, 227], [89, 159, 218],
                      [0, 83, 170], [0, 4, 167], [0, 140, 1]]) / 255.0

# A. Standard PR experiments
for val, col in zip(PR_values, pr_colors):
    if val in [0, 10, 30]:
        base_dir = f"/data92/PeterChang/back_to_master1220/Moist_Dycore/IdealizeSpetral.jl/exp/HSt42/HSt42_{val}/"
    elif val in [20, 40, 50]:
        base_dir = f"/data92/PeterChang2/Moist_Dycore/IdealizeSpetral.jl/exp/HSt42/HSt42_{val}/"
        
    experiments.append({
        "label": f"PR{val}",
        "L_val": val * 0.01,
        "color": col,
        "linestyle": "-",
        "raw_base_dir": base_dir,
        "file_prefix": f"RH80_PR{val}_20000day_startfrom_"
    })

# B. Midlat experiment
experiments.append({
    "label": "Midlat",
    "L_val": 0.5, 
    "color": "red",
    "linestyle": "--",
    "raw_base_dir": "/data92/PeterChang/Moist_Dycore_only_midlatitude_LH/IdealizeSpetral.jl/exp/HSt42/HSt42_50/",
    "file_prefix": "RH80_PR50_20000day_startfrom_" 
})

# C. Trop experiment
experiments.append({
    "label": "Trop",
    "L_val": 0.5, 
    "color": "orange",
    "linestyle": "--",
    "raw_base_dir": "/data92/PeterChang2/Moist_Dycore_only_tropical_LH_30S_30N/IdealizeSpetral.jl/exp/HSt42/HSt42_50/",
    "file_prefix": "RH80_PR50_20000day_startfrom_" 
})

# ==========================================
# 2. Helper Functions
# ==========================================
def get_lanczos_weights(cutoff_period_days=10, window_days=20):
    """Generates 41 Lanczos weights for daily data."""
    fc = 1.0 / cutoff_period_days 
    k = np.arange(-window_days, window_days + 1)
    weights = 2 * fc * np.sinc(2 * fc * k) * np.sinc(k / window_days)
    return weights / weights.sum()

def compute_convergence(flux, yd, cy, radius=EARTH_RADIUS):
    """
    Calculates spherical convergence using explicit finite differences.
    flux shape: (Time, Z, Lat)
    """
    dmdy = np.zeros_like(flux)
    num_lat = flux.shape[2]
    
    # Central difference for interior points
    for k in range(1, num_lat - 1):
        num = (flux[:, :, k+1] * (cy[k+1]**2)) - (flux[:, :, k-1] * (cy[k-1]**2))
        den = radius * (cy[k]**2) * (yd[k+1] - yd[k-1])
        dmdy[:, :, k] = -num / den
        
    # Equator Boundary (index 0)
    # Safe because cy[0] is cos(0) = 1.0
    dmdy[:, :, 0] = -((flux[:, :, 1] * (cy[1]**2)) - (flux[:, :, 0] * (cy[0]**2))) / \
                     (radius * (cy[0]**2) * (yd[1] - yd[0]))
                     
    # Pole Boundary (index -1)
    # cy[-1] is cos(90) = 0.0. Division by zero causes a singularity.
    # Physically, momentum flux vanishes at the pole, so we set convergence to 0.
    dmdy[:, :, -1] = 0.0
                      
    return dmdy

# Generate latitude arrays (Assuming 0 to 90 N for 32 points)
lat_deg = np.linspace(0, 90, Y_DIM)
lat_rad = np.deg2rad(lat_deg)
cy = np.cos(lat_rad)
area_weights = cy ** 0.5

# Filter and Spectral parameters
weights = get_lanczos_weights(cutoff_period_days=10, window_days=30)
nperseg = 256  # 256 day sections 
noverlap = 128 # 128 day overlap 
fs = 1.0       # 1 cycle/day

# IMPORTANT: Make sure your custom EOF module is imported here if it's an external script!
# e.g., from my_eof_module import EOF

# ==========================================
# 3. Main Loop Over Experiments
# ==========================================
for exp_idx, exp in enumerate(experiments):
    label = exp["label"]
    print(f"\n{'='*40}")
    print(f"Processing Experiment {exp_idx + 1}/{len(experiments)}: {label}")
    print(f"{'='*40}")
    
    # 3.1 Load Data and Compute Daily Means
    u_prime_full = np.zeros((total_days, Z_DIM, Y_DIM, X_DIM), dtype=np.float32)
    v_prime_full = np.zeros((total_days, Z_DIM, Y_DIM, X_DIM), dtype=np.float32)
    
    # Array to store the daily, zonal-mean zonal wind [u] for EOF calculation
    Uzm = np.zeros((total_days, Z_DIM, Y_DIM), dtype=np.float32)
    
    current_idx = 0
    for day in start_days:
        # Check extensions just in case (e.g., .dat vs .h5)
        filename = os.path.join(exp["raw_base_dir"], f"{exp['file_prefix']}{day}day_final.dat") 
        if not os.path.exists(filename):
            filename = os.path.join(exp["raw_base_dir"], f"{exp['file_prefix']}{day}day_final.h5")
        
        with h5py.File(filename, 'r') as f:
            raw_u_6h = f["grid_u_c_xyzt"][:, :, -Y_DIM:, :]
            raw_v_6h = f["grid_v_c_xyzt"][:, :, -Y_DIM:, :]
            
            raw_u_daily = raw_u_6h.reshape(DAYS_PER_FILE, 4, Z_DIM, Y_DIM, X_DIM).mean(axis=1)
            raw_v_daily = raw_v_6h.reshape(DAYS_PER_FILE, 4, Z_DIM, Y_DIM, X_DIM).mean(axis=1)
            
            # Calculate zonal means
            u_zonal_mean = raw_u_daily.mean(axis=-1)
            v_zonal_mean = raw_v_daily.mean(axis=-1)
            
            # Store [u] for EOF calculation
            Uzm[current_idx : current_idx + DAYS_PER_FILE] = u_zonal_mean
            
            # Calculate eddy component: u' = u - [u]
            u_prime_full[current_idx : current_idx + DAYS_PER_FILE] = raw_u_daily - np.expand_dims(u_zonal_mean, axis=-1)
            v_prime_full[current_idx : current_idx + DAYS_PER_FILE] = raw_v_daily - np.expand_dims(v_zonal_mean, axis=-1)
            
        current_idx += DAYS_PER_FILE
        if day % 2500 == 0:
            print(f"  Processed up to day {day}...")

    print(f"[{label}] Data loaded. Applying Lanczos filter...")
    
    # # 3.2 Time Filtering
    # u_prime_l = convolve1d(u_prime_full, weights, axis=0, mode='reflect')
    # v_prime_l = convolve1d(v_prime_full, weights, axis=0, mode='reflect')
    # u_prime_h = u_prime_full - u_prime_l
    # v_prime_h = v_prime_full - v_prime_l
    
    # # FREE MEMORY IMMEDIATELY
    # del u_prime_full, v_prime_full 

    # # 3.3 Flux and Convergence
    # print(f"[{label}] Computing momentum fluxes and convergence...")
    # synoptic_flux = (u_prime_h * v_prime_h).mean(axis=-1)
    # residual_flux = (u_prime_l * v_prime_l).mean(axis=-1)
    
    # conv_synoptic = compute_convergence(synoptic_flux, lat_rad, cy)
    # conv_residual = compute_convergence(residual_flux, lat_rad, cy)
    # 3.2 Time Filtering
    u_prime_l = convolve1d(u_prime_full, weights, axis=0, mode='reflect')
    v_prime_l = convolve1d(v_prime_full, weights, axis=0, mode='reflect')
    
    u_prime_h = u_prime_full - u_prime_l
    v_prime_h = v_prime_full - v_prime_l
    
    # NEW: Calculate the TOTAL eddy flux BEFORE deleting the full arrays
    # Taking the zonal mean immediately so it uses almost no memory
    total_flux = (u_prime_full * v_prime_full).mean(axis=-1)
    
    # FREE MEMORY IMMEDIATELY
    del u_prime_full, v_prime_full 

    # 3.3 Flux and Convergence
    print(f"[{label}] Computing momentum fluxes and convergence...")
    
    # 1. Synoptic Flux [u'_h * v'_h]
    synoptic_flux = (u_prime_h * v_prime_h).mean(axis=-1)
    
    # 2. Residual Flux is Total minus Synoptic (per LH2001 definition)
    residual_flux = total_flux - synoptic_flux
    
    # Calculate convergence (dmdy) as before
    conv_synoptic = compute_convergence(synoptic_flux, lat_rad, cy)
    conv_residual = compute_convergence(residual_flux, lat_rad, cy)
    
    # 3.4 Vertical Average and Anomalies
    m_synoptic_2d = conv_synoptic.mean(axis=1)  
    m_residual_2d = conv_residual.mean(axis=1)  

    m_synoptic_anom = m_synoptic_2d - m_synoptic_2d.mean(axis=0)
    m_residual_anom = m_residual_2d - m_residual_2d.mean(axis=0)

    # ==========================================
    # 3.5 Calculate EOF for THIS Experiment
    # ==========================================
    print(f"[{label}] Projecting onto EOF1...")
    
    # Take the vertical average of the zonal-mean wind to get <[u]> (Time, Y)
    Uzm_vert_avg = Uzm.mean(axis=1) 
    
    # Compute anomalies by subtracting the time-mean
    # Uzm_vert_avg_anom = Uzm_vert_avg - Uzm_vert_avg.mean(axis=0)
    
    # Weight by the square root of cosine latitude
    u_weighted = Uzm_vert_avg * (cy ** 0.5) 
    
    # Run the EOF calculation
    single_EOF_u = EOF((u_weighted,), n_components=Y_DIM, field="1D")
    single_EOF_u.get()
    
    # Extract the leading EOF spatial pattern
    EOF_u_final = single_EOF_u.EOF[0]
    PC1 = single_EOF_u.PC[0]
    
    EOF_u_final_normalized = EOF_u_final / EOF_u_final.std()
    
    # 3.6 Projection
    m_synoptic_t = np.dot(m_synoptic_anom * area_weights, EOF_u_final_normalized) 
    m_residual_t = np.dot(m_residual_anom * area_weights, EOF_u_final_normalized)

    # 3.7 Spectral Analysis
    print(f"[{label}] Computing power spectra...")
    freq_syn, power_syn = welch(m_synoptic_t, fs=fs, window='hann', nperseg=nperseg, noverlap=noverlap)
    freq_res, power_res = welch(m_residual_t, fs=fs, window='hann', nperseg=nperseg, noverlap=noverlap)

    # ==========================================
    # 3.8 Save Results to HDF5
    # ==========================================
    save_path = f"/data92/PeterChang/back_to_master1220/Moist_Dycore/IdealizeSpetral.jl/exp/HSt42/eddy_fluxes_and_spectra_{label}.h5"
    
    print(f"[{label}] Saving data to {save_path}...")
    with h5py.File(save_path, 'w') as f:
        # 1. Save 1D Power Spectra arrays
        f.create_dataset("freq_syn", data=freq_syn)
        f.create_dataset("P_synoptic", data=power_syn)
        f.create_dataset("freq_res", data=freq_res)
        f.create_dataset("P_residual", data=power_res)
        
        # 2. Save the 1D m(t) time series
        f.create_dataset("m_synoptic_t", data=m_synoptic_t)
        f.create_dataset("m_residual_t", data=m_residual_t)
        
        # 3. Save the full 3D zonal-mean momentum fluxes [u'v'] (Time, Z, Y)
        f.create_dataset("synoptic_flux", data=synoptic_flux, compression="gzip", compression_opts=4)
        f.create_dataset("residual_flux", data=residual_flux, compression="gzip", compression_opts=4)

        # 4. Save the full 3D eddy momentum flux convergence (dmdy) (Time, Z, Y)
        f.create_dataset("conv_synoptic", data=conv_synoptic, compression="gzip", compression_opts=4)
        f.create_dataset("conv_residual", data=conv_residual, compression="gzip", compression_opts=4)

        # Optional: Save the EOF pattern and Uzm anomaly if you want to reuse them later
        f.create_dataset("EOF1_pattern", data=EOF_u_final_normalized)
        f.create_dataset("PC1", data=PC1)

    print(f"[{label}] Saved successfully!")

print("\nAll 8 experiments processed successfully!")