* NODESET is a temporary nonlinear constraint and must be released before convergence.
.OPTIONS NODESET_ITERS=3 NODESET_G=1e5
VDD vdd 0 DC 1.8
R1 vdd out 10k
D1 out 0 DTEST
.MODEL DTEST D IS=1e-14 N=1
.NODESET V(out)=0.65
.OP
.END
