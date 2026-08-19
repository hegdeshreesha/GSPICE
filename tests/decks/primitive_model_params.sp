* Primitive diode/BJT fidelity parameters accepted by parser and solver
Vd in 0 DC 1
Rd in dout 10
D1 dout 0 DMOD
.MODEL DMOD D IS=1e-14 N=1 RS=25 BV=0.5 IBV=1u NBV=1.2 CJO=0.1p

Vc c 0 DC 1
Vb b 0 DC 0.72
Q1 c b 0 QMOD
.MODEL QMOD NPN IS=1e-16 BF=100 BR=2 NF=1 NR=1 IKF=1m IKR=1m VAF=50 VAR=20 ISE=1e-13 ISC=1e-13 NE=1.5 NC=2 CJE=0.1p CJC=0.05p TF=1p

.TRAN 1n 1n
.SAVE I(D1) I(Q1)
.MEASURE TRAN id FIND I(D1) AT=1n
.MEASURE TRAN iq FIND I(Q1) AT=1n
.END
