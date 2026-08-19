* Semiconductor current probe convention: BJT current is collector terminal current.
.MODEL QN NPN IS=1e-16 BF=100 NF=1
VCC vcc 0 DC 5
VB base 0 DC 0.7
RC vcc out 1k
Q1 out base 0 QN
.TRAN 1n 1n
.SAVE I(Q1)
.MEAS TRAN iq1 FIND I(Q1) AT=1n
.END
