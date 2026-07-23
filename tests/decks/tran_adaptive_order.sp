* Smooth variable-step BDF run should exercise p-1/p/p+1 order control.
.OPTIONS METHOD=BDF MAXORD=5 ADAPTIVE=1 LTE_MODE=PREDICTOR LTE_AUDIT_INTERVAL=8 RELTOL=1e-6 VNTOL=1e-10
V1 in 0 DC 1
R1 in out 1k
C1 out 0 1n
.TRAN 1n 100n 0 10n
.END
