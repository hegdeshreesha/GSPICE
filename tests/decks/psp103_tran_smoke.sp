* PSP103 transient + AC smoke through GSDI/GMC
.MODEL nch psp103va TYPE=1 VTO=0.45 BET=1m LAMBDA=0.02
VDD vdd 0 DC 1.8
VIN gate 0 PULSE(1.2 0.3 0 1n 1n 10n 20n)
RD vdd drain 10k
CL drain 0 10f
M1 drain gate 0 0 nch W=1u L=0.13u
.TRAN 1n 60n
.PRINT TRAN v(drain)
.END
