* Common PSP aliases should be consumed, while truly unsupported params remain visible.
.MODEL nch psp103va TYPE=1 VFB0=0.45 PHIBO=0.7 BETA=1m LAMDA=0.02 TOX=2n
VDD vdd 0 DC 1.8
VIN gate 0 DC 1.2
RD vdd drain 10k
M1 drain gate 0 0 nch W=1u L=0.13u
.OP
.END
