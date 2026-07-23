* State-preserving OSDI F/Q/Jacobian and charge-conservation audit.
.PRE_OSDI "../../osdi/psp103.osdi"
.MODEL pch psp103va type=-1
.OPTIONS DAE_AUDIT=1 DAE_AUDIT_TOL=2e-3
VDD vdd 0 DC 2
VG gate 0 DC 0
RLOAD out 0 100k
N1 out gate vdd vdd pch w=1u l=0.13u nf=1 mult=1
.OP
.END
