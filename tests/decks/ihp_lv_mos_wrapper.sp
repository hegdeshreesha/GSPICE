* IHP low-voltage MOS wrapper fallback smoke deck
.SUBCKT sg13_lv_nmos d g s b w=0.15u l=0.13u ng=1 m=1
.ENDS
.SUBCKT sg13_lv_pmos d g s b w=0.15u l=0.13u ng=1 m=1
.ENDS
XM1 out in vdd vdd sg13_lv_pmos l=0.13u w=0.15u ng=1 m=1
XM2 out in 0 0 sg13_lv_nmos l=0.13u w=0.15u ng=1 m=1
VDD vdd 0 DC 2
VIN in 0 DC 0 PULSE(0 2 0 1n 1n 1u 2u)
.OPTIONS ADAPTIVE=0
.OP
.TRAN 1n 10n
.END
