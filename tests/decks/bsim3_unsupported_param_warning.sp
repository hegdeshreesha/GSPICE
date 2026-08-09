* Unsupported BSIM3 model parameters are verbose-only diagnostics.
.MODEL nref nmos level=49 vth0=0.4 toxe=10n PHANTOMVAR=1.0
VDS d 0 0.8
VGS g 0 0.7
M1 d g 0 0 nref W=1u L=1u
.OP
.END
