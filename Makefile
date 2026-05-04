# Usage:
#   make LANG=C [target]       - build C version
#   make LANG=Fortran [target] - build Fortran version (default)
#
# Targets: all, lib, example, clean

LANG ?= Fortran

ifeq ($(LANG),C)
    SUBDIR = 00_C
else ifeq ($(LANG),Fortran)
    SUBDIR = 01_Fortran
else
    $(error Unsupported LANG=$(LANG). Use LANG=C or LANG=Fortran)
endif

TARGETS = all lib example clean

.PHONY: $(TARGETS)

$(TARGETS):
	$(MAKE) -C $(SUBDIR) $@
