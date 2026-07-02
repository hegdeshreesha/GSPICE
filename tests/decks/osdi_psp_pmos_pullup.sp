* Real OpenVAF PSP PMOS pull-up regression.
.PRE_OSDI "../../osdi/psp103.osdi"
.MODEL pch psp103va type=-1
VDD vdd 0 DC 2
VG gate 0 DC 0
RLOAD out 0 100k
N1 out gate vdd vdd pch w=1u l=0.13u nf=1 mult=1
.OP
.END
