* .OPTIONS SAVE=SELECTED should write only explicit .SAVE signals.
.OPTIONS SAVE=SELECTED
V1 in 0 PULSE(0 1 0 100p 100p 1n 2n)
R1 in out 1k
C1 out 0 1p
.SAVE V(out)
.TRAN 1n 2n
.END
