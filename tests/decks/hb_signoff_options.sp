* Global signoff mode should require native HB and tighten iteration limits.
.OPTIONS SIGNOFF=YES
V1 in 0 SIN(0 1 1k)
R1 in out 1k
C1 out 0 1n
.HB 1k NHARMS=1 MAXITER=12 TOL=100
.END
