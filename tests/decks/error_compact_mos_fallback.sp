Compact MOS fallback must be explicit
.MODEL sg13_lv_nmos psp103va(type=1)
M1 out in 0 0 sg13_lv_nmos L=0.13u W=0.15u
VDD out 0 1.0
VIN in 0 1.0
.OP
.END
