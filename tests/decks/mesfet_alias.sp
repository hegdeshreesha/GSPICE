MESFET alias routed through primitive JFET approximation
VDD vdd 0 DC 5
VGG gate 0 DC -1
RD vdd drain 1k
.MODEL NM NMES(BETA=1m VTO=-2 LAMBDA=0.02 IS=1e-14)
J1 drain gate 0 NM
.END
