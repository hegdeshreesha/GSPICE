* Semiconductor current probe convention: diode current is anode-to-cathode.
.MODEL DMOD D IS=1e-14 N=1
V1 an 0 DC 0.7
D1 an 0 DMOD
.TRAN 1n 1n
.SAVE I(D1)
.MEAS TRAN id1 FIND I(D1) AT=1n
.END
