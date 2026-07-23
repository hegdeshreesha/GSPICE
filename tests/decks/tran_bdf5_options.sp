* General variable-step BDF integration smoke deck
.OPTIONS METHOD=BDF MAXORD=5 ADAPTIVE=0 RELTOL=1e-6 VNTOL=1e-10
V1 in 0 PULSE(0 1 0 5p 5p 100p 250p)
R1 in out 1k
C1 out 0 20f
.TRAN 10p 200p
.END
