Parameter expression regression
.param base = 1k
.param scale=2
.param rval={base*scale}
V1 in 0 DC {2*scale}
R1 in out {rval}
R2 out 0 {base}
.OP
.END
