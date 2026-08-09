* IHP SG13G2 Schottky diode wrapper through primitive diode/BJT/passive support.
.LIB "C:/EDA/LumenCircuitStudio/external/ihp_pdk/ihp-sg13g2/libs.tech/ngspice/models/cornerDIO.lib" dio_tt
VA a 0 DC 0.4
VC c 0 DC 0
XS1 a c 0 schottky_nbl1 l=1u w=0.3u Nx=1 Ny=1
.OP
.END
