* Native HB PSP103 smoke: driven gate produces periodic drain response.
.MODEL nch psp103va TYPE=1 VTO=0.45 BET=1m LAMBDA=0.02
VDD vdd 0 DC 1.8
VIN gate 0 SIN(0.9 0.05 1k)
RD vdd drain 10k
CL drain 0 10f
M1 drain gate 0 0 nch W=1u L=0.13u
.HB 1k NHARMS=2 MAXITER=12 TOL=100 SIGNOFF=YES
.END
