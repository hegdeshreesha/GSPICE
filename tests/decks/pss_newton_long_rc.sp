* PSS finite-difference shooting Newton smoke deck.
V1 in 0 SIN(0 1 1k)
R1 in out 10k
C1 out 0 1u
.OPTIONS RELTOL=1e-6 VNTOL=1u ABSTOL=1p
.SAVE V(out)
.PSS 1k 1 TSTAB=1m PSS_ADAPTIVE=YES MAX_PSS_ITER=5 PSS_RESIDUAL_GOAL=1e3
.END
