* PMOS LEVEL=103 implies PSP polarity even when TYPE is omitted.
.MODEL pch PMOS LEVEL=103 VTO=0.45 BET=1m LAMBDA=0.02
VDD vdd 0 DC 1.8
VIN gate 0 DC 0.6
RD drain 0 10k
M1 drain gate vdd vdd pch W=1u L=0.13u
.OP
.END
