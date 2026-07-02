* Lumen/IHP wrapper smoke deck for manual PSP verification.
.LIB "C:\EDA\LumenCircuitStudio\external\ihp_pdk\ihp-sg13g2\libs.tech\ngspice\models\cornerMOSlv.lib" mos_tt
XM1 out in vdd vdd sg13_lv_pmos l=0.13u w=0.15u ng=1 m=1
XM3 out in 0 0 sg13_lv_nmos l=0.13u w=0.15u ng=1 m=1
V1 vdd 0 DC 2
V3 in 0 DC 0 PULSE(0 2 0 1n 1n 1u 2u)
.OP
.TRAN 100n 2u 0
.END
