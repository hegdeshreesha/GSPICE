* Native GSDI/GMC PSP103 first-light operating-point deck.
.GSDI "psp103_native.gsdi"
.MODEL nch psp103va TYPE=1 VTO=0.45 BET=1m LAMBDA=0.02
VDD vdd 0 DC 1.8
VIN gate 0 DC 1.2
RD vdd drain 10k
N1 drain gate 0 0 nch W=1u L=0.13u
.OP
.END
