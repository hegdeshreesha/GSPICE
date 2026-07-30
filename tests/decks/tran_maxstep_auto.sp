* MAXSTEP=AUTO should ignore the SPICE .TRAN tmax field.
.OPTIONS MAXSTEP=AUTO
V1 in 0 PULSE(0 1 0 100p 100p 1n 2n)
R1 in out 1k
C1 out 0 1p
.TRAN 1n 2n 0 20p
.END
