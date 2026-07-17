Lumen IHP inverter transient short regression
.LIB "C:\EDA\LumenCircuitStudio\external\ihp_pdk\ihp-sg13g2\libs.tech\ngspice\models\cornerMOSlv.lib" mos_tt
XM1 net2 net0 net1 net1 sg13_lv_pmos l=0.13u w=0.15u ng=1 m=1
XM3 net2 net0 0 0 sg13_lv_nmos l=0.13u w=0.15u ng=1 m=1
V1 net1 0 DC 2
V3 net0 0 DC 0 PULSE(0 2 0 10n 10n 1u 2u)
C0 net2 0 100f
.OPTIONS ACCURACY=HIGH METHOD=AUTO ADAPTIVE=1 RELTOL=3e-4 VNTOL=300n ABSTOL=100f TRTOL=1e-3 TRABSTOL=300n ITL4=80
.OP
.TRAN 50p 50n 0 5p
.END
