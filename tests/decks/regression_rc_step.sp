RC step numeric regression
.OPTIONS ADAPTIVE=0 METHOD=BE RELTOL=1e-5 VNTOL=1e-9
V1 in 0 DC 0 PULSE(0 1 0 1n 1n 2m 4m)
R1 in out 1k
C1 out 0 1u
.TRAN 10u 1m 0
.END
