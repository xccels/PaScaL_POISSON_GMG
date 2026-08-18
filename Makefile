# Usage:
#   make LANG=C [target]       - build C version
#   make LANG=Fortran [target] - build Fortran version (default)
#
# Targets
#   LANG=C       : all (= cpu), cpu, gpu, clean, clean_cpu, clean_gpu
#   LANG=Fortran : all, lib, example, clean

# LANG is also the POSIX locale environment variable, which make imports as a
# make variable. Ignore that inherited value so a bare `make` still defaults to
# Fortran; an explicit `make LANG=C` on the command line still wins.
ifeq ($(origin LANG),environment)
    LANG = Fortran
endif
LANG ?= Fortran

ifeq ($(LANG),C)
    SUBDIR = 00_C
else ifeq ($(LANG),Fortran)
    SUBDIR = 01_Fortran
else
    $(error Unsupported LANG=$(LANG). Use LANG=C or LANG=Fortran)
endif

TARGETS = all lib example clean cpu gpu clean_cpu clean_gpu

.PHONY: $(TARGETS)

$(TARGETS):
	$(MAKE) -C $(SUBDIR) $@
