* Runtime DAE derivative and charge-conservation audit
.OPTIONS DAE_AUDIT=1 DAE_AUDIT_TOL=1e-3
VDD d 0 1.2
VG g 0 1.0
R1 d out 1k
C1 out 0 10f
D1 out 0 DTEST
M1 out g 0 0 NM W=2u L=1u
.MODEL DTEST D(IS=1e-14 N=1 CJO=2f)
.MODEL NM NMOS(VTO=0.45 KP=120u)
.OP
.END
