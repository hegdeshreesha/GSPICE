BJT primitive regression
.MODEL QN NPN(IS=1e-16 BF=100 NF=1)
VCC vcc 0 DC 5
VB base 0 DC 0.7
RC vcc out 1k
Q1 out base 0 QN
.OP
.END
