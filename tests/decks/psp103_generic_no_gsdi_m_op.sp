* Generic PSP103 MOS syntax should not need an IHP wrapper or explicit .GSDI.
.MODEL nch psp103va TYPE=1 VTO=0.45 BET=1m LAMBDA=0.02
VDD vdd 0 DC 1.8
VIN gate 0 DC 1.2
RD vdd drain 10k
M1 drain gate 0 0 nch W=1u L=0.13u
.OP
.END
