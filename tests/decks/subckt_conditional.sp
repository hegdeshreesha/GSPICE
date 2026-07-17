Subcircuit conditional regression
.subckt gated in out sel=1 rval=1k
.if (sel == 1)
R1 in out {rval}
.else
R1 in out {rval*2}
.endif
.ends
V1 in 0 1
X1 in out gated sel=1 rval=1k
Rload out 0 1k
.OP
.END
