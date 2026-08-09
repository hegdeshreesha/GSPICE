* PSP NQS family spelling should route without an explicit .GSDI.
.MODEL nch pspnqs103e TYPE=1 VTO=0.45 BET=1m LAMBDA=0.02
VDD vdd 0 DC 1.8
VIN gate 0 DC 1.2
RD vdd drain 10k
N1 drain gate 0 0 nch W=1u L=0.13u
.OP
.END
