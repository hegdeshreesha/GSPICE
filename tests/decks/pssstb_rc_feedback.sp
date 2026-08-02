* PSSSTB two-port smoke deck.
V1 in 0 DC 1
R1 in loop_in 1k
WPROBE loop_in loop_out
EGAIN out 0 loop_out 0 -10
RLOAD out 0 10k
RFB out loop_out 10k
CFB loop_out 0 1n
.PSSSTB 1k DEC 3 10 10k SIDEBANDS=1
.END
