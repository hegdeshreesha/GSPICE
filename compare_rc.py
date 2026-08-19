import re, sys, math

def parse_gspice_csv(path):
    rows = {}
    with open(path) as f:
        header = f.readline().strip().split(',')
        name = header[1].strip('"') if len(header) > 1 else 'v'
        for line in f:
            p = line.strip().split(',')
            if len(p) >= 2:
                rows[float(p[0])] = float(p[1])
    return rows, name

def parse_ngspice(path):
    rows = {}
    with open(path, errors='replace') as f:
        for line in f:
            parts = line.split('\t')
            if len(parts) >= 3:
                try:
                    t = float(parts[1]); v = float(parts[2])
                    rows[t] = v
                except ValueError:
                    pass
    return rows, 'v(out)'

def parse_xyce_prn(path):
    rows = {}
    with open(path) as f:
        header = f.readline()
        cols = header.strip().split()
        idx_v = [i for i, c in enumerate(cols) if c.upper().startswith('V(')][0]
        for line in f:
            p = line.split()
            if len(p) > idx_v:
                try:
                    t = float(p[1]); rows[t] = float(p[idx_v])
                except ValueError:
                    pass
    return rows, 'V(out)'

def compare(a, b, label, tol=1e-6):
    keys = sorted(set(a) & set(b))
    max_rel, max_abs, worst_t = 0.0, 0.0, None
    for t in keys:
        d = abs(a[t] - b[t])
        denom = max(abs(a[t]), abs(b[t]), 1e-12)
        r = d / denom
        if r > max_rel: max_rel, worst_t = r, t
        if d > max_abs: max_abs = d
    print(f"{label:45s} n={len(keys):5d}  max_rel={max_rel:.3e} max_abs={max_abs:.3e} worst_t={worst_t}")
    return max_rel

def compare_interp(a, b, label, exclude=()):
    ts = sorted(t for t in a if not any(lo <= t <= hi for lo, hi in exclude))
    bt = sorted(b.keys())
    max_rel, max_abs, worst_t = 0.0, 0.0, None
    n = 0
    j = 0
    for t in ts:
        while j < len(bt) - 2 and bt[j + 1] <= t:
            j += 1
        t0, t1 = bt[j], bt[j + 1]
        if t1 == t0:
            continue
        v = b[t0] + (b[t1] - b[t0]) * (t - t0) / (t1 - t0)
        d = abs(a[t] - v)
        denom = max(abs(a[t]), abs(v), 1e-12)
        r = d / denom
        n += 1
        if r > max_rel: max_rel, worst_t = r, t
        if d > max_abs: max_abs = d
    print(f"{label:45s} n={n:5d}  max_rel={max_rel:.3e} max_abs={max_abs:.3e} worst_t={worst_t}")
    return max_rel

g, gn = parse_gspice_csv('rc_gspice.raw')
n, nn = parse_ngspice('rc_ngspice.out')
x, xn = parse_xyce_prn('rc_x.prn')
print(f"points: gspice={len(g)} ngspice={len(n)} xyce={len(x)}")
r1 = compare_interp(g, n, 'GSPICE vs ngspice (TRAN RC, excl 1ms edge)', exclude=((0.00095, 0.00106),))
r2 = compare_interp(g, x, 'GSPICE vs Xyce    (TRAN RC, excl 1ms edge)', exclude=((0.00095, 0.00106),))
r3 = compare_interp(n, x, 'ngspice vs Xyce   (TRAN RC, excl 1ms edge)', exclude=((0.00095, 0.00106),))