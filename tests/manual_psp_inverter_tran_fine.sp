* Fine-step PSP inverter transient diagnostic.
.PRE_OSDI "C:\EDA\GSPICE\osdi\psp103.osdi"
.MODEL nch psp103va type=1
.MODEL pch psp103va type=-1
VDD vdd 0 DC 2
VIN in 0 DC 0 PULSE(0 2 0 1n 1n 5n 10n)
NP out in vdd vdd pch w=1u l=0.13u nf=1 mult=1
NN out in 0 0 nch w=1u l=0.13u nf=1 mult=1
.TRAN 0.2n 20n 0
.END
