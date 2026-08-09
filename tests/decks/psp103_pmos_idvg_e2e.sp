* PSP103 native PMOS via M-prefix device
.option temp=21
.include psp103_pmos-2.mod
vd d 0 dc -0.1
vg g 0 dc 0.0
vs s 0 dc 0.0
vb b 0 dc 0.0
mp1 d g s b pch l=0.1u w=1u m=1
.dc Vg 0 -1.5 -0.05
.print dc v(g) i(vd)
.end