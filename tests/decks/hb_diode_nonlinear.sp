* Native HB nonlinear smoke: diode conduction produces a visible second harmonic.
V1 in 0 SIN(0.55 0.08 1k)
R1 in out 100
D1 out 0 DHB
.MODEL DHB D(IS=1e-14 N=1 CJO=0.1p)
.HB 1k NHARMS=2 MAXITER=12 TOL=100 SIGNOFF=YES
.END
