* Two instances share one OSDI .MODEL block but retain private instance state.
.PRE_OSDI "../../osdi/psp103.osdi"
.MODEL pch psp103va type=-1
.OPTIONS DAE_AUDIT=1 DAE_AUDIT_TOL=2e-3
VDD vdd 0 DC 2
VG gate 0 DC 0
R1 out1 0 100k
R2 out2 0 200k
N1 out1 gate vdd vdd pch w=1u l=0.13u nf=1 mult=1
N2 out2 gate vdd vdd pch w=2u l=0.13u nf=1 mult=1
.OP
.END
