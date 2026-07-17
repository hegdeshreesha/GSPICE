IHP primitive fallback must be explicit
.SUBCKT sg13_lv_nmos d g s b w=0.15u l=0.13u ng=1 m=1
.ENDS
XM1 out in 0 0 sg13_lv_nmos l=0.13u w=0.15u ng=1 m=1
VIN in 0 DC 1
.OP
.END
