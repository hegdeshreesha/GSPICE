* Minimal OpenVAF PSP smoke deck for manual verification.
.PRE_OSDI "C:\EDA\GSPICE\osdi\psp103.osdi"
.MODEL nch psp103va type=1
VDD drain 0 DC 1.0
VIN gate 0 DC 0.8
N1 drain gate 0 0 nch w=1u l=0.13u nf=1 mult=1
.OP
.END
