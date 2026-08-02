* Explicit current saves are written alongside save-all voltages.
V1 in 0 DC 1
R1 in out 1k
C1 out 0 1p
.TRAN 1n 2n
.SAVE ALL
.SAVE I(V1)
.SAVE I(R1)
.END
