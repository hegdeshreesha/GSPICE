* OSDI PSP transient-noise colored-source smoke deck.
.PRE_OSDI "../../osdi/psp103.osdi"
.MODEL pch psp103va type=-1
VDD vdd 0 DC 2
VG gate 0 DC 0
RLOAD out 0 100k
N1 out gate vdd vdd pch w=1u l=0.13u nf=1 mult=1
.OPTIONS TRAN_PROGRESS_INTERVAL=0
.TRANNOISE 1n 3n FMIN=1 FMAX=1e6 SEED=77 SCALE=0.001 NOISEMODE=ZOH TONES_PER_DEC=2
.END
