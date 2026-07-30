* Nonlinear diode LPTV PNOISE smoke deck.
V1 in 0 SIN(0 0.45 1k)
R1 in out 200
D1 out 0 DMIX
C1 out 0 10n
.MODEL DMIX D(IS=1e-14 N=1.2 CJO=0.5p)
.PNOISE V(out) V1 DEC 2 100 1k FUND=1k SIDEBANDS=1
.END
