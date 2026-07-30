* PSS driven shooting residual smoke deck.
V1 in 0 SIN(0 1 1k)
R1 in out 1k
C1 out 0 1n
.OPTIONS RELTOL=1e-6 VNTOL=1u ABSTOL=1p
.SAVE V(out)
.PSS 1k 1 TSTAB=1m PSS_ADAPTIVE=YES MAX_PSS_ITER=10
.END
