import time
import numpy as np

# Simulate 500,000 point waveform dataset
N = 500000
x_arr = np.linspace(0, 1.1e-6, N)
y_arr = np.sin(2 * np.pi * 57.5e6 * x_arr) + 0.1 * np.random.randn(N)

x_min, x_max = 0.0, 1.1e-6
y_min, y_max = -1.5, 1.5
width_px = 1200
height_px = 600
left_px = 78
bottom_px = 548

def ultra_vectorized_decimate_and_map(x_arr, y_arr, x_min, x_max, y_min, y_max, left_px, bottom_px, width_px, height_px):
    n = len(x_arr)
    if n == 0:
        return np.array([])
    
    # 1. Searchsorted slice
    i0 = max(0, np.searchsorted(x_arr, x_min, side='left') - 1)
    i1 = min(n, np.searchsorted(x_arr, x_max, side='right') + 1)
    
    n_vis = i1 - i0
    target_bins = width_px * 2
    
    if n_vis <= target_bins:
        x_sub = x_arr[i0:i1]
        y_sub = y_arr[i0:i1]
    else:
        bin_len = n_vis // target_bins
        if bin_len < 1:
            x_sub = x_arr[i0:i1]
            y_sub = y_arr[i0:i1]
        else:
            tot = target_bins * bin_len
            x_g = x_arr[i0:i0+tot].reshape(target_bins, bin_len)
            y_g = y_arr[i0:i0+tot].reshape(target_bins, bin_len)
            
            min_i = np.argmin(y_g, axis=1)
            max_i = np.argmax(y_g, axis=1)
            rows = np.arange(target_bins)
            
            idx1 = np.minimum(min_i, max_i)
            idx2 = np.maximum(min_i, max_i)
            
            x_sub = np.empty(target_bins * 2, dtype=np.float64)
            y_sub = np.empty(target_bins * 2, dtype=np.float64)
            
            x_sub[0::2] = x_g[rows, idx1]
            x_sub[1::2] = x_g[rows, idx2]
            y_sub[0::2] = y_g[rows, idx1]
            y_sub[1::2] = y_g[rows, idx2]

    # 2. Vectorized Screen Coordinate Transformation (SIMD)
    dx = (x_max - x_min) or 1.0
    dy = (y_max - y_min) or 1.0
    
    sx = left_px + ((x_sub - x_min) / dx) * width_px
    sy = bottom_px - ((y_sub - y_min) / dy) * height_px
    
    # Interleave into (N, 2) screen points array
    pts = np.column_stack((sx, sy))
    return pts

# Benchmark 500,000 points
t0 = time.perf_counter()
for _ in range(1000):
    pts = ultra_vectorized_decimate_and_map(x_arr, y_arr, x_min, x_max, y_min, y_max, left_px, bottom_px, width_px, height_px)
t1 = time.perf_counter()

print(f"Ultra-vectorized time per frame (500,000 points): {(t1 - t0)/1000*1000:.4f} ms")
print(f"Output screen points: {len(pts)}")
