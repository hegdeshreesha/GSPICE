import time
from bisect import bisect_left, bisect_right
import numpy as np

# Simulate a 100,000 point waveform dataset
N = 100000
x_data = np.linspace(0, 1.1e-6, N).tolist()
y_data = np.sin(2 * np.pi * 57.5e6 * np.array(x_data)).tolist()

def fast_viewport_render_samples(x_data: list[float], y_data: list[float], x_min: float, x_max: float, width_px: int) -> tuple[list[float], list[float]]:
    n = min(len(x_data), len(y_data))
    if n == 0:
        return [], []
    
    # 1. O(log N) Viewport Slicing
    i_start = bisect_left(x_data, x_min)
    if i_start > 0:
        i_start -= 1
    i_end = bisect_right(x_data, x_max)
    if i_end < n:
        i_end += 1
        
    n_vis = i_end - i_start
    target_bins = max(300, int(width_px) * 2)
    
    # If visible points are fewer than target_bins, no decimation needed!
    if n_vis <= target_bins:
        return x_data[i_start:i_end], y_data[i_start:i_end]
        
    # 2. Fast MinMax Pyramidal Binning (M-4 algorithm)
    bin_size = n_vis / target_bins
    out_x = []
    out_y = []
    
    for b in range(target_bins):
        b0 = i_start + int(b * bin_size)
        b1 = i_start + int((b + 1) * bin_size)
        if b0 >= b1:
            continue
        
        # Sliced range
        sub_y = y_data[b0:b1]
        
        # Find min & max in this pixel bin
        min_idx = 0
        max_idx = 0
        min_val = sub_y[0]
        max_val = sub_y[0]
        
        for k in range(1, len(sub_y)):
            v = sub_y[k]
            if v < min_val:
                min_val = v
                min_idx = k
            if v > max_val:
                max_val = v
                max_idx = k
                
        # Append in chronological order
        idx1, idx2 = (min_idx, max_idx) if min_idx <= max_idx else (max_idx, min_idx)
        
        out_x.append(x_data[b0 + idx1])
        out_y.append(y_data[b0 + idx1])
        if idx1 != idx2:
            out_x.append(x_data[b0 + idx2])
            out_y.append(y_data[b0 + idx2])
            
    return out_x, out_y

# Test zoom into 1% window vs full window
t0 = time.perf_counter()
for _ in range(500):
    # Zoomed into 1% window (100ns to 110ns)
    rx, ry = fast_viewport_render_samples(x_data, y_data, 100e-9, 110e-9, 1200)
t1 = time.perf_counter()

print(f"Zoomed 1% window render time per frame: {(t1 - t0)/500*1000:.4f} ms (Points: {len(rx)})")

t0 = time.perf_counter()
for _ in range(500):
    # Full window (0 to 1.1us)
    rx, ry = fast_viewport_render_samples(x_data, y_data, 0, 1.1e-6, 1200)
t1 = time.perf_counter()

print(f"Full window render time per frame: {(t1 - t0)/500*1000:.4f} ms (Points: {len(rx)})")
