Transient measure INTEG and DERIV regression
.OPTIONS ADAPTIVE=0 METHOD=BE RELTOL=1e-9 VNTOL=1e-12
V1 out 0 PWL(0 0 1m 1)
R1 out 0 1k
.TRAN 10u 1m 0
.MEAS TRAN area INTEG V(out) FROM=0 TO=1m
.MEAS TRAN slope DERIV V(out) FROM=0 TO=1m
.MEAS TRAN slope_at DERIV V(out) AT=0.5m
.END
