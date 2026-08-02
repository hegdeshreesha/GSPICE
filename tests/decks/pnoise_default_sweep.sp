* PNOISE should default the offset sweep when only an output node is provided.
V1 in 0 SIN(0 1 1k)
R1 in out 1k
C1 out 0 1n
.PSS 1k 1
.PNOISE V(out)
.END
