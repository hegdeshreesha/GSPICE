GSPICE PSP103.4 inverter transient reference
.PRE_OSDI "../../osdi/psp103.osdi"
.MODEL nch psp103va type=1
.MODEL pch psp103va type=-1
.OPTIONS ACCURACY=HIGH ADAPTIVE=1 METHOD=GEAR2 LTE_RELTOL=1e-3 TRABSTOL=300n TRTOL=1 CHGTOL=1e-15 MAXORD=2
VDD vdd 0 DC 2
VIN in 0 DC 0 PULSE(0 2 100p 20p 20p 200p 500p)
NP out in vdd vdd pch w=1u l=0.13u nf=1 mult=1
NN out in 0 0 nch w=1u l=0.13u nf=1 mult=1
CLOAD out 0 20f
.TRAN 2p 1n 0 5p
.END
