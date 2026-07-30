High accuracy transient preset regression
.OPTIONS ACCURACY=HIGH TRAN_PROGRESS_INTERVAL=0
V1 in 0 DC 0 PULSE(0 1 0 1n 1n 5n 10n)
R1 in out 1k
C1 out 0 1p
.TRAN 1n 10n
.END
