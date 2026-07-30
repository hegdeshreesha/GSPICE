* PNOISE shifted-sideband smoke deck.
V1 in 0 DC 0 AC 1
R1 in out 1k
C1 out 0 1n
.PNOISE V(out) V1 DEC 3 10 10k FUND=1k SIDEBANDS=1
.END
