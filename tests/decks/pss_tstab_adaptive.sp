* PSS tstab and adaptive-orbit smoke deck.
V1 in 0 SIN(0 1 1k)
R1 in out 1k
C1 out 0 1u
.SAVE V(out)
.PSS 1k 1 TSTAB=1m PSS_ADAPTIVE=YES PSS_CONTINUATION=NO
.END
