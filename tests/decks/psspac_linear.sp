* PSSPAC shifted-sideband small-signal smoke deck.
V1 in 0 DC 0 AC 1
R1 in out 1k
C1 out 0 1n
.SAVE V(out)
.PSSPAC 1k DEC 3 10 10k SIDEBANDS=1
.END
