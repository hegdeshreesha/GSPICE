* BSIM4 state reference for reduced-model diagnostics.
.model NREF nmos level=54 vth0=0.40 u0=500 toxe=10n k1=0 k2=0 dvt0=0 dvt1=0.53 dvt2=-0.032
+ xj=0.15u ua=2.0e-9 ub=5.0e-19 vsat=1.0e5
VDS d 0 0.80
VGS g 0 0
M1 d g 0 0 NREF W=1u L=1u
.dc VGS 0 1 0.1
.print dc v(g) @m1[vth] @m1[vdsat] @m1[gm] @m1[gds]
.end
