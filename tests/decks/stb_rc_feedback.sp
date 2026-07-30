* First-generation STB return-ratio smoke deck.
* WPROBE is a zero-volt loop-break probe in DC and a 1 V return-ratio
* excitation during .STB.
V1 in 0 DC 1
R1 in loop_in 1k
WPROBE loop_in loop_out
EGAIN out 0 loop_out 0 -10
RLOAD out 0 10k
RFB out loop_out 10k
CFB loop_out 0 1n
.STB DEC 5 1 1e6
.END
