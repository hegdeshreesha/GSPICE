* Transient method/options parser and engine smoke deck
.OPTIONS METHOD=GEAR2 MAXSTEP=20p ACCURACY=HIGH
V1 in 0 PULSE(0 1 0 100p 100p 1n 2n)
R1 in out 1k
C1 out 0 20f
.TRAN 200p 2n
.END
