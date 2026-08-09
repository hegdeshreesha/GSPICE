* Native JUNCAP Express parser/DAE smoke test
.model JX JUNCAPEXP SWJUNEXP=1 PHITD=0.026 ISATFOR1=1e-14 MFOR1=1 CJO=1e-12
V1 in 0 PULSE(0 0.6 0 1n 1n 1n 4n)
D1 in 0 JX
.op
.tran 1n 4n
.end
