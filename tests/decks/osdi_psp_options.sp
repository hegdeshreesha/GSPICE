* OSDI advanced options should be available through SPICE .OPTIONS
.PRE_OSDI "../../osdi/psp103.osdi"
.MODEL nch psp103va TYPE=1
.MODEL pch psp103va TYPE=-1
VDD vdd 0 DC 2
VIN in 0 PULSE(0 2 0 10n 10n 100n 200n)
NP out in vdd vdd pch w=0.15u l=0.13u
NN out in 0 0 nch w=0.15u l=0.13u
C0 out 0 10f
.OPTIONS OSDI_LIMITING_RHS=1 OSDI_SPICE_RHS=1 OSDI_BIND_FULL_MODEL_PARAMS=0 ADAPTIVE=0 METHOD=BE
.TRAN 1n 20n 0 1n
.SAVE V(out) V(in) V(vdd)
.END
