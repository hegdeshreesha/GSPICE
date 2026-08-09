* IHP SG13G2 LV PSP AC smoke.
.LIB "C:/EDA/LumenCircuitStudio/external/ihp_pdk/ihp-sg13g2/libs.tech/ngspice/models/cornerMOSlv.lib" mos_tt
VDD vdd 0 DC 1.8
VG gate 0 DC 1.0 AC 1
RD vdd drain 10k
CL drain 0 1f
XM1 drain gate 0 0 sg13_lv_nmos l=0.13u w=0.15u ng=1 m=1
.AC LIN 5 1k 1MEG
.END
