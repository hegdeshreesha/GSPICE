DC sweep and primitive model-card regression
.MODEL DNOM D(IS=1e-14 N=1.0 CJO=2p)
V1 in 0 DC 0
R1 in out 1k
D1 out 0 DNOM
.DC V1 0 1 0.5
.END
