* IHP SG13G2 high-voltage MOS varicap AC C-V smoke.
.LIB "C:/EDA/LumenCircuitStudio/external/ihp_pdk/ihp-sg13g2/libs.tech/ngspice/models/cornerMOShv.lib" mos_tt
VG1 g1 0 DC -1 AC 1
VG2 g2 0 DC 0
XCV g1 well g2 0 sg13_hv_svaricap l=600n w=3u Nx=2 Ny=1
.AC LIN 4 1k 1MEG
.END
