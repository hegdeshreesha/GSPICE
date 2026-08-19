import re

def parse_ngspice_raw(path):
    vals = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(('Title', 'Date', 'Plotname', 'Flags', 'No.', 'Variables', 'Values')) or line == '0':
                continue
            parts = line.split()
            try:
                data = [float(p) for p in parts]
            except ValueError:
                continue
            if len(data) == 8:
                vals.append(data)
    return vals

def parse_gspice_dc_stdout(path):
    rows = {}
    with open(path, encoding='utf-16', errors='ignore') as f:
        for line in f:
            m = re.match(r'\s*([\d.eE+-]+)\s*\|\s*([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+(-?[\d.eE+-]+)\s*$', line)
            if m:
                vg = float(m.group(1))
                ivd = float(m.group(6))
                rows[vg] = ivd
    return rows

raw = parse_ngspice_raw(r"C:\EDA\_psp_ref_work\psp_ref_idvg.raw")
csv = parse_gspice_dc_stdout(r"C:\EDA\GSPICE\psp_idvg_gspice_stdout.txt")
print(f"ngspice rows={len(raw)} gspice rows={len(csv)}")
worst = (0, None, None)
rows = []
for r in raw:
    vg = r[0]
    if vg in csv:
        id_ref = -r[7]
        id_gs = csv[vg]
        rel = abs(id_gs - id_ref) / max(abs(id_ref), 1e-15)
        absd = abs(id_gs - id_ref)
        rows.append((vg, id_gs, id_ref, rel, absd))
        if abs(id_ref) > 1e-12 and rel > worst[0]:
            worst = (rel, vg, (id_gs, id_ref))
print(f"{'Vg':>8} {'Id_gspice':>16} {'Id_ngspiceOSDI':>16} {'rel':>10} {'abs':>10}")
for vg, a, b, rel, abd in rows:
    print(f"{vg:8.3f} {a:16.6e} {b:16.6e} {rel:10.3e} {abd:10.3e}")
print(f"\nWORST (rel, |Id|>1e-12): {worst}")