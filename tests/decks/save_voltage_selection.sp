* .SAVE voltage selection should restrict transient RAW variables.
V1 in 0 DC 0 PULSE(0 1 0 1n 1n 5n 10n)
R1 in out 1k
C1 out 0 1p
.TRAN 1n 20n 0 1n
.SAVE V(out)
.END
