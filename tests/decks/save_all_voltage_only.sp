* Save-all writes node voltages, not branch currents.
V1 in 0 DC 1
R1 in out 1k
C1 out 0 1p
.TRAN 1n 2n
.SAVE ALL
.END
