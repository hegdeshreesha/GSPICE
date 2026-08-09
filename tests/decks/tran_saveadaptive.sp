* SAVEADAPTIVE should preserve accepted internal transient points in RAW.
V1 in 0 PULSE(0 1 0 100p 100p 5n 10n)
R1 in out 1k
C1 out 0 1p
.OPTIONS RELTOL=1e-4 VNTOL=1e-7 ABSTOL=1e-12 MAXSTEP=100p SAVEADAPTIVE=1
.TRAN 1n 20n
.SAVE V(out)
.END
