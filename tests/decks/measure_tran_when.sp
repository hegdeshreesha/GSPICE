Transient measure WHEN crossing regression
.OPTIONS ADAPTIVE=0 METHOD=BE RELTOL=1e-9 VNTOL=1e-12
V1 out 0 PWL(0 0 1m 1)
R1 out 0 1k
.TRAN 10u 1m 0
.MEAS TRAN t50 WHEN V(out)=0.5 RISE=1
.END
