* OSDI parser/loader/device plumbing smoke deck
.OSDI builtin:mos_level_50
.MODEL nch mos_level_50 type=1
VDD vdd 0 DC 1.8
VIN gate 0 DC 0.8
RD vdd drain 10k
N1 drain gate 0 0 nch w=1u l=0.13u
.OP
.END
