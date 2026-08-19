* Native HBNOISE nonlinear smoke: signoff HB orbit feeds diode periodic noise conversion.
V1 in 0 SIN(0.55 0.08 1k)
R1 in out 100
D1 out 0 DHB
.MODEL DHB D(IS=1e-14 N=1 CJO=0.1p)
.HBNOISE V(out) 1k LIN 1 1k 1k SIDEBANDS=2 SIGNOFF=YES
.END
