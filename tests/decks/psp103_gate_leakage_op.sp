* PSP gate leakage knobs should contribute native gate current.
.MODEL nch psp103va TYPE=1 VTO=0.45 BET=1m IGINV=1u IGOV=0.5u IGOVD=0.25u
VDD vdd 0 DC 1.8
VIN gate 0 DC 1.2
RD vdd drain 10k
M1 drain gate 0 0 nch W=2u L=0.13u
.OP
.END
