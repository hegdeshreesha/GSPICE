* Scanner regression: accepted-but-approximate PSP103 parameters must be enforceable.
.MODEL nch psp103va TYPE=1 VTO=0.45 BET=1m SWNUD=1
VDD d 0 DC 1
VG g 0 DC 1
M1 d g 0 0 nch W=1u L=0.13u
.OP
.END
