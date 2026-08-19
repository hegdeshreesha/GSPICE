* .FOUR transient Fourier post-analysis smoke
V1 out 0 SIN(0 1 1k)
R1 out 0 1k
.TRAN 10u 5m
.FOUR 1k V(out)
.END
