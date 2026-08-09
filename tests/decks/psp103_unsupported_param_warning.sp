* Unsupported PSP parameters should be visible, not silently swallowed.
.MODEL nch psp103va TYPE=1 VTO=0.45 BET=1m RSH=10 XJ=20n TOX=2n
VDD vdd 0 DC 1.8
VIN gate 0 DC 1.2
RD vdd drain 10k
M1 drain gate 0 0 nch W=1u L=0.13u
.OP
.END
