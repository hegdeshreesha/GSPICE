* PSP103 demo: foundry model file (CEA 103.3 ref data), bias sweep,
* tempo corner, and a transient stepped gate. All through GSDI + GMC.
.INCLUDE psp103_nmos-2.mod
VDD vdd 0 DC 1.8
RB vdd d 10k
M1 d g 0 0 nch L=0.1u W=1u M=4
VG g 0 DC 0.9
.DC VG 0 1.2 0.05
.TRAN 1n 10n
.AC DEC 10 1k 10g