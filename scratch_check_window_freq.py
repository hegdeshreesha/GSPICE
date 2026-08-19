import numpy as np
from scratch_check_freq import read_raw

var_names, data = read_raw("test_ring.raw")
t = data[:, 0]
stg3 = data[:, var_names.index("V(STG3)")]
outnet = data[:, var_names.index("V(Outnet)")]

# 1. Measurement over early startup window (0 to 100ns) vs full window
def window_freq(t, v, t_min, t_max):
    mask = (t >= t_min) & (t <= t_max)
    t_sub = t[mask]
    v_sub = v[mask]
    crossings = []
    for i in range(len(v_sub)-1):
        if v_sub[i] <= 0.5 and v_sub[i+1] > 0.5:
            tc = t_sub[i] + (0.5 - v_sub[i]) * (t_sub[i+1] - t_sub[i]) / (v_sub[i+1] - v_sub[i])
            crossings.append(tc)
    if len(crossings) >= 2:
        return 1.0 / np.mean(np.diff(crossings)) / 1e6
    return 0.0

print("\n--- Frequency vs Time Window ---")
print(f"STG3   (0 - 50ns startup) : {window_freq(t, stg3, 0, 50e-9):.2f} MHz")
print(f"Outnet (0 - 50ns startup) : {window_freq(t, outnet, 0, 50e-9):.2f} MHz")
print(f"STG3   (100ns - 1.1us steady): {window_freq(t, stg3, 100e-9, 1.1e-6):.2f} MHz")
print(f"Outnet (100ns - 1.1us steady): {window_freq(t, outnet, 100e-9, 1.1e-6):.2f} MHz")

# 2. Discrete index distance without sub-sample interpolation
def discrete_index_freq(v, dt=22e-12):
    peaks = []
    for i in range(1, len(v)-1):
        if v[i-1] < v[i] and v[i] >= v[i+1] and v[i] > 0.8:
            peaks.append(i)
    if len(peaks) >= 2:
        index_diffs = np.diff(peaks)
        avg_indices = np.mean(index_diffs)
        return 1.0 / (avg_indices * dt) / 1e6
    return 0.0

print("\n--- Frequency via Discrete Peak Sample Index Counting ---")
print(f"STG3   (discrete sample steps): {discrete_index_freq(stg3):.2f} MHz")
print(f"Outnet (discrete sample steps): {discrete_index_freq(outnet):.2f} MHz")
