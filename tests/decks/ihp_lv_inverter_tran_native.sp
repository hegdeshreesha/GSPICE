* IHP SG13G2 LV PSP transient smoke.
.LIB "C:/EDA/LumenCircuitStudio/external/ihp_pdk/ihp-sg13g2/libs.tech/ngspice/models/cornerMOSlv.lib" mos_tt
VDD vdd 0 DC 1.8
VIN in 0 DC 0 PULSE(0 1.8 0 20p 20p 200p 500p)
XM1 out in vdd vdd sg13_lv_pmos l=0.13u w=0.3u ng=1 m=1
XM2 out in 0 0 sg13_lv_nmos l=0.13u w=0.15u ng=1 m=1
CL out 0 10f
.OPTIONS METHOD=AUTO ADAPTIVE=1 RELTOL=3e-4 VNTOL=300n ABSTOL=100f TRTOL=1 LTE_RELTOL=1e-3 TRABSTOL=300n ITL4=80
.TRAN 20p 2n 0 5p
.END
