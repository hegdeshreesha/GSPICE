Initial condition, temperature, and global compatibility regression
.GLOBAL vdd vss
.TEMP 85
VDD vdd 0 DC 1.8
R1 vdd out 1k
C1 out 0 1p
.IC V(out)=0.7
.TRAN 1n 2n 0 1n UIC
.END
