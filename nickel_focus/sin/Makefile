ifndef LROOT
	LROOT = /usr/local/lick
endif

include $(LROOT)/etc/Make.defs

SINFILES = $(wildcard *.sin)
FILES = $(SINFILES:%.sin=%)

BINFILES = $(FILES)
BINSUB = nickel

include $(LROOT)/etc/Make.rules
