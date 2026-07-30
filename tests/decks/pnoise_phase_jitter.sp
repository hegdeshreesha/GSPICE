* PNOISE phase-noise and jitter flow smoke deck.
V1 in 0 SIN(0 1 1k)
R1 in out 1k
C1 out 0 1u
.PNOISE V(out) V1 DEC 2 100 1k FUND=1k SIDEBANDS=1 PHASENOISE=YES JITTER=YES CARRIER=1k
.END
