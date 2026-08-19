* Semiconductor current probe convention: MOS current is drain-to-source.
.MODEL nch NMOS VTO=0.5 KP=1m LAMBDA=0.02
VD d 0 DC 1.8
VG g 0 DC 1.2
M1 d g 0 0 nch W=1u L=1u
.TRAN 1n 1n
.SAVE I(M1)
.MEAS TRAN im1 FIND I(M1) AT=1n
.END
