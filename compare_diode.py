import re

def parse_ng(path):
    rows = {}
    with open(path, errors='replace') as f:
        for line in f:
            parts = line.split('\t')
            if len(parts) >= 3:
                try:
                    rows[float(parts[1])] = float(parts[2])
                except ValueError:
                    pass
    return rows

def parse_xyce(path):
    rows = {}
    with open(path) as f:
        f.readline()
        for i, line in enumerate(f):
            p = line.split()
            if len(p) >= 2:
                try:
                    rows[round(i * 0.01, 6)] = float(p[-1])
                except ValueError:
                    pass
    return rows

def parse_gspice(path):
    rows = {}
    with open(path, encoding='utf-16', errors='ignore') as f:
        for line in f:
            m = re.match(r'\s*([\d.eE+-]+)\s*\|\s*([\d.eE+-]+)\s+([\d.eE+-]+)\s*$', line)
            if m:
                rows[float(m.group(1))] = abs(float(m.group(3)))
    return rows

n = parse_ng('diode_ng.out')
x = parse_xyce('diode_x.prn')
g = parse_gspice('diode_gspice_stdout.txt')
print(f"ngspice={len(n)} xyce={len(x)} gspice={len(g)}")
for key, name in ((n,'ngspice'),(x,'Xyce'),(g,'GSPICE')):
    v = key.get(0.3)
    print(f"Id(0.3V) {name}: {v:.6e}")
print("\nVd        gspice          ngspice         xyce           g/n      g/x")
for vd in [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:
    gn, nn, xn = g.get(vd), n.get(vd), x.get(vd)
    if gn and nn and xn:
        print(f"{vd:.1f}  {gn:12.6e}  {nn:12.6e}  {xn:12.6e}  {abs(gn/nn-1):8.2e}  {abs(gn/xn-1):8.2e}")