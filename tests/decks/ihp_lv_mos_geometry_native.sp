* IHP SG13G2 LV PSP geometry sweep smoke.
.LIB "C:/EDA/LumenCircuitStudio/external/ihp_pdk/ihp-sg13g2/libs.tech/ngspice/models/cornerMOSlv.lib" mos_tt
VDD vdd 0 DC 1.8
VG gate 0 DC 1.0
RD1 vdd d1 10k
RD2 vdd d2 10k
RD3 vdd d3 10k
XM1 d1 gate 0 0 sg13_lv_nmos l=0.13u w=0.15u ng=1 m=1
XM2 d2 gate 0 0 sg13_lv_nmos l=0.26u w=0.3u ng=2 m=1
XM3 d3 gate 0 0 sg13_lv_nmos l=1u w=1u ng=4 m=2
.OP
.END
