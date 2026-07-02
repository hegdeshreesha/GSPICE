* Subckt and simple global parameter smoke deck
.PARAM RV=1k
.SUBCKT rcell a b
R1 a b RV
.ENDS
V1 in 0 1
X1 in out rcell
R2 out 0 1k
.OP
.TRAN 1n 2n
.END
