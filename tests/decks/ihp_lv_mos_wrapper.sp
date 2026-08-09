* IHP low-voltage MOS wrapper native PSP smoke deck
.GSDI "psp103_native.gsdi"
.SUBCKT sg13_lv_nmos d g s b w=0.15u l=0.13u ng=1 m=1
.ENDS
.SUBCKT sg13_lv_pmos d g s b w=0.15u l=0.13u ng=1 m=1
.ENDS
.MODEL sg13g2_lv_nmos_psp psp103va TYPE=1 VTO=0.45 BET=1m LAMBDA=0.02
VDD vdd 0 DC 1.8
VIN in 0 DC 1.2
RD vdd out 10k
XM2 out in 0 0 sg13_lv_nmos l=0.13u w=0.15u ng=1 m=1
.OP
.END
