Monte Carlo production spec regression
V1 in 0 DC 1
R1 in out 1k
R2 out 0 1k
.SPEC vout V(out) MIN=0.49 MAX=0.51
.MC 4 V1 GAUSS(1 0) SEED=11
.END
