* PSS initial-condition option parser smoke deck.
V1 in 0 SIN(0 1 1k)
R1 in out 1k
C1 out 0 1u
.IC out=0.25
.SAVE V(out)
.PSS 1k 1 USE_INITIAL_CONDITIONS=YES
.END
