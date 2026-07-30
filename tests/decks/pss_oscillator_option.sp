* PSS oscillator option parser smoke deck.
V1 in 0 SIN(0 1 1k)
R1 in out 1k
C1 out 0 1u
.SAVE V(out)
.PSS 1k 1 OSCILLATOR=YES PSS_ADAPTIVE=YES MAX_PSS_ITER=12 PSS_RESIDUAL_GOAL=10
.END
