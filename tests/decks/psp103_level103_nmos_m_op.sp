* PSP decks may spell PSP103 as NMOS LEVEL=103.
.MODEL nch NMOS LEVEL=103 VTO=0.45 BET=1m LAMBDA=0.02
VDD vdd 0 DC 1.8
VIN gate 0 DC 1.2
RD vdd drain 10k
M1 drain gate 0 0 nch W=1u L=0.13u
.OP
.END
