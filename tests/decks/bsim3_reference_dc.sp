* Independent BSIM3 reference deck for ngspice comparison.
* Level 49 is the ngspice BSIM3 model family.
.model NREF nmos level=49 vth0=0.40 kp=120u gamma=0.50 phi=0.60
+ u0=500 tox=10n xj=0.15u
VDS d 0 0.80
VGS g 0 0
M1 d g 0 0 NREF W=1u L=1u
.dc VGS 0 1 0.1
.print dc v(g) i(VDS)
.end
