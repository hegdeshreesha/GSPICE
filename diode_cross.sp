* Diode DC sweep cross-check
vd d 0 dc 0
.model diodep D(IS=1e-14 N=1.5 RS=0 CJO=0)
d1 d 0 diodep
.dc vd 0 0.9 0.01
.print dc i(vd)
.end