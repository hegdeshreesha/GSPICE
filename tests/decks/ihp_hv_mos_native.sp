* IHP SG13G2 high-voltage MOS wrapper through native PSP path.
.LIB "C:/EDA/LumenCircuitStudio/external/ihp_pdk/ihp-sg13g2/libs.tech/ngspice/models/cornerMOShv.lib" mos_tt
VDD vdd 0 DC 3.3
VIN gate 0 DC 2.0
RD vdd drain 10k
XM1 drain gate 0 0 sg13_hv_nmos l=0.34u w=0.35u ng=1 m=1
.OP
.END
