* PSP GIDL knobs should contribute drain/bulk leakage.
.MODEL nch psp103va TYPE=1 VTO=0.45 BET=1m AGIDLD=1m BGIDLD=0.1
VDD vdd 0 DC 1.8
VIN gate 0 DC 0
RD vdd drain 10k
M1 drain gate 0 0 nch W=2u L=0.13u
.OP
.END
