* Cross-simulator diode DC reference deck.
.MODEL DMOD D IS=1e-14 N=1
VSWEEP in 0 DC 0
R1 in out 1k
D1 out 0 DMOD
.DC VSWEEP 0 0.8 0.2
.PRINT DC V(out)
.END
