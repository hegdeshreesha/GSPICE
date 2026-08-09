* IHP SG13G2 MIM capacitor wrapper with model TC1/TC2 temperature scaling.
.LIB "C:/EDA/LumenCircuitStudio/external/ihp_pdk/ihp-sg13g2/libs.tech/ngspice/models/cornerCAP.lib" cap_typ
.TEMP 125
V1 in 0 DC 1
XC1 in 0 cap_cmim l=7u w=7u
.OP
.END
