* Unsupported BSIM4 model parameters must be visible, not silently swallowed.
.MODEL nref nmos level=54 vth0=0.4 toxe=10n w0=0.5u PHANTOMVAR=1.0
VDS d 0 0.8
VGS g 0 0.7
M1 d g 0 0 nref W=1u L=1u
.OP
.END