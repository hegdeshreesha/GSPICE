* BSIM3 must not silently become a primitive MOS model.
.model NREF nmos level=49 vth0=0.40 kp=120u
VDS d 0 0.8
VGS g 0 0.8
M1 d g 0 0 NREF W=1u L=1u
.op
.end
