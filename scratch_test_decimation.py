import time
import numpy as np

n_pts = 100000
x_arr = np.linspace(0, 1.1e-6, n_pts)
y_arr = np.sin(2 * np.pi * 57.5e6 * x_arr)

x_min, x_max = 100e-9, 200e-9
width_px = 1200

def vectorized_viewport_minmax(x_data, y_data, x_min, x_max, width_px=1200):
    # 1. Binary search slice for visible viewport
    i_start = max(0, np.searchsorted(x_data, x_min, side='left') - 1)
    i_end = min(len(x_data), np.searchsorted(x_data, x_max, side='right') + 1)
    
    count = i_end - i_start
    target_bins = width_px * 2
    if count <= target_bins:
        return x_data[i_start:i_end], y_data[i_start:i_end]
    
    # 2. Reshape into bins and compute min/max using vector operations
    x_sub = x_data[i_start:i_end]
    y_sub = y_data[i_start:i_end]
    
    # Trim to multiple of target_bins
    bin_len = count // target_bins
    if bin_len < 1:
        return x_sub, y_sub
    
    trimmed_len = target_bins * bin_len
    x_grid = x_sub[:trimmed_len].reshape(target_bins, bin_len)
    y_grid = y_sub[:trimmed_len].reshape(target_bins, bin_len)
    
    min_idx = np.argmin(y_grid, axis=1)
    max_idx = np.argmax(y_grid, axis=1)
    
    rows = np.arange(target_bins)
    
    # Pair indices in chronological order
    idx1 = np.minimum(min_idx, max_idx)
    idx2 = np.maximum(min_idx, max_idx)
    
    res_x = np.empty(target_bins * 2)
    res_y = np.empty(target_bins * 2)
    
    res_x[0::2] = x_grid[rows, idx1]
    res_x[1::2] = x_grid[rows, idx2]
    res_y[0::2] = y_grid[rows, idx1]
    res_y[1::2] = y_grid[rows, idx2]
    
    return res_x, res_y

t0 = time.perf_counter()
for _ in range(1000):
    rx, ry = vectorized_viewport_minmax(x_arr, y_arr, x_min, x_max, width_px)
t1 = time.perf_counter()

print(f"Vectorized Viewport MinMax time per frame: {(t1 - t0) / 1000 * 1000:.4f} ms")
print(f"Decimated output points: {len(rx)}")
