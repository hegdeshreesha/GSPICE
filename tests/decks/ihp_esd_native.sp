* IHP SG13G2 ESD diode wrapper through primitive diode support.
.LIB "C:/EDA/LumenCircuitStudio/external/ihp_pdk/ihp-sg13g2/libs.tech/ngspice/models/cornerDIO.lib" dio_tt
VDD vdd 0 DC 1.8
VPAD pad 0 DC 0.9
XD1 vdd pad 0 diodevdd_2kv m=1
XD2 vdd pad 0 diodevss_2kv m=1
.OP
.END
