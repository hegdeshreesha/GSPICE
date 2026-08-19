* Independent BSIM3 reference deck for ngspice comparison.
* Level 49 is the ngspice BSIM3 model family.
.model NREF nmos level=49 vth0=0.40 u0=500 tox=10n xj=0.15u
+ vsat=1e5 nfactor=1.2 k1=0 k2=0 ua=2.0e-9 ub=5.0e-19
VDS d 0 0.80
VGS g 0 0.2
M1 d g 0 0 NREF W=1u L=1u
.dc VGS 0.2 1.0 0.1
.print dc v(g) i(VDS)
.end
