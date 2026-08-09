* IHP SG13G2 low-voltage MOS library first-light native PSP deck.
.LIB "C:/EDA/LumenCircuitStudio/external/ihp_pdk/ihp-sg13g2/libs.tech/ngspice/models/cornerMOSlv.lib" mos_tt
VDD vdd 0 DC 1.8
VIN gate 0 DC 1.2
RD vdd drain 10k
XM1 drain gate 0 0 sg13_lv_nmos l=0.13u w=0.15u ng=1 m=1
.OP
.END
