Transient measure TRIG/TARG delay regression
.OPTIONS ADAPTIVE=0 METHOD=BE RELTOL=1e-9 VNTOL=1e-12
VIN in 0 PWL(0 0 1m 1)
VOUT out 0 PWL(0 0 0.2m 0 1.2m 1)
RIN in 0 1k
ROUT out 0 1k
.TRAN 10u 1.2m 0
.MEAS TRAN delay TRIG V(in) VAL=0.5 RISE=1 TARG V(out) VAL=0.5 RISE=1
.END
