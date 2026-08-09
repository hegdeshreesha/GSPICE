* Native JUNCAP2 TAT + avalanche parser/GSDI smoke test.
.model J2 JUNCAP2 PHITD=0.026 PHITR=0.026 IDSAT=1e-14 CSRH=1e-13 CTAT=1e-18 CBBT=1e-18 VBR=0.2 XJUN=1u CJORBOT=1e-3
V1 in 0 -0.5
D1 in 0 J2
.op
.end
