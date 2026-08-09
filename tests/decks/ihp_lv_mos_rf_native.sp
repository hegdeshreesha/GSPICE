* IHP SG13G2 LV PSP RF model-card smoke.
.LIB "C:/EDA/LumenCircuitStudio/external/ihp_pdk/ihp-sg13g2/libs.tech/ngspice/models/cornerMOSlv.lib" mos_tt
VDD vdd 0 DC 1.8
VG gate 0 DC 1.0 AC 1
RD vdd drain 10k
M1 drain gate 0 0 sg13g2_lv_nmos_psp_rf l=0.13u w=0.15u ng=1 m=1
.OP
.AC LIN 3 1k 1MEG
.END
