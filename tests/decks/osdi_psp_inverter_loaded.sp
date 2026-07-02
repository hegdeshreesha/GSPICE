* Real OpenVAF PSP inverter transient regression with output load.
.PRE_OSDI "../../osdi/psp103.osdi"
.MODEL nch psp103va type=1
.MODEL pch psp103va type=-1
VDD vdd 0 DC 2
VIN in 0 DC 0 PULSE(0 2 0 20p 20p 100p 240p)
NP out in vdd vdd pch w=1u l=0.13u nf=1 mult=1
NN out in 0 0 nch w=1u l=0.13u nf=1 mult=1
CLOAD out 0 20e-15
.TRAN 5p 300p 0
.END
