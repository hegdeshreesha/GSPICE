* PSP loaded-gate DC operating point regression.
* A low input must pull the first inverter output high even when it drives
* another PSP gate. This requires expanded hidden OSDI internal equations.
.PRE_OSDI "../../osdi/psp103.osdi"
.MODEL nch psp103va type=1
.MODEL pch psp103va type=-1
VDD vdd 0 DC 2
VIN in 0 DC 0
NP1 out in vdd vdd pch w=1u l=0.13u nf=1 mult=1
NN1 out in 0 0 nch w=1u l=0.13u nf=1 mult=1
NP2 load out vdd vdd pch w=1u l=0.13u nf=1 mult=1
NN2 load out 0 0 nch w=1u l=0.13u nf=1 mult=1
.OP
.END
