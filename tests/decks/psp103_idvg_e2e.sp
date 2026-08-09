* PSP103 native IdVg end-to-end parity (ngspice reference: ref_idvg.sp, temp=21)
.option temp=21
.include psp103_nmos-2.mod
vd d 0 dc 0.1
vg g 0 dc 0.0
vs s 0 dc 0.0
vb b 0 dc 0.0
nm1 d g s b nch l=0.1u w=1u m=1
.dc Vg 0 1.5 0.05
.print dc v(g) i(vd)
.end