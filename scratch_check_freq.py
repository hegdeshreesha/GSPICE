import numpy as np

def read_raw(filename):
    with open(filename, 'r', encoding='latin-1') as f:
        lines = f.readlines()
    
    num_vars = 0
    num_pts = 0
    var_names = []
    values_start = 0
    
    for i, line in enumerate(lines):
        line_str = line.strip()
        if line_str.startswith("No. Variables:"):
            num_vars = int(line_str.split(":")[1].strip())
        elif line_str.startswith("No. Points:"):
            num_pts = int(line_str.split(":")[1].strip())
        elif line_str.startswith("Variables:"):
            for j in range(i + 1, i + 1 + num_vars):
                parts = lines[j].strip().split()
                if len(parts) >= 2:
                    var_names.append(parts[1])
        elif line_str.startswith("Values:"):
            values_start = i + 1
            break
            
    print(f"Num Vars: {num_vars}, Num Points: {num_pts}")
    print(f"Var Names: {var_names}")
    
    data = []
    cur_pt = []
    for line in lines[values_start:]:
        parts = line.strip().split()
        for p in parts:
            try:
                cur_pt.append(float(p))
                if len(cur_pt) == num_vars:
                    data.append(cur_pt)
                    cur_pt = []
            except ValueError:
                continue
                
    data = np.array(data)
    return var_names, data

var_names, data = read_raw("test_ring.raw")
t = data[:, 0]

def calc_freq(t, v, thresh=0.5):
    crossings = []
    for i in range(len(v)-1):
        if v[i] <= thresh and v[i+1] > thresh:
            tc = t[i] + (thresh - v[i]) * (t[i+1] - t[i]) / (v[i+1] - v[i])
            crossings.append(tc)
    
    if len(crossings) < 2:
        return 0.0, 0, []
    
    periods = np.diff(crossings)
    steady_periods = periods[int(len(periods)*0.2):]
    avg_period = np.mean(steady_periods)
    freq = 1.0 / avg_period
    return freq, len(crossings), steady_periods

print("\n--- Measured Frequencies (50% VDD threshold) ---")
for idx, name in enumerate(var_names):
    if idx == 0: continue
    freq, count, periods = calc_freq(t, data[:, idx], 0.5)
    std_p = np.std(periods) if len(periods) > 0 else 0
    print(f"{name:<12}: Freq = {freq/1e6:8.4f} MHz | Count = {count:4d} | Period StdDev = {std_p*1e12:6.2f} ps")

