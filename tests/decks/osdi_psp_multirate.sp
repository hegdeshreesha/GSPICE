* OSDI PSP transient must converge when MULTIRATE is enabled.
.PRE_OSDI "../../osdi/psp103.osdi"
.MODEL nch psp103va type=1
.MODEL pch psp103va type=-1
.OPTIONS ADAPTIVE=1 MULTIRATE=1 METHOD=GEAR2 LTE_RELTOL=1e-3 TRABSTOL=300n TRTOL=1
VDD vdd 0 DC 1.2
VIN in 0 PULSE(0 1.2 0 10p 10p 200p 400p)
NP out in vdd vdd pch w=2u l=0.045u
NN out in 0 0 nch w=1u l=0.045u
.TRAN 10p 400p
.END
