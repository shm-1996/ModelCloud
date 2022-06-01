#
#  Makefile for Semenov et al. (2003) opacities
#

FCOMP = mpif90 -fdefault-real-8 -fdefault-double-8
# FCOMP = openmpif90 -fdefault-real-8 -fdefault-double-8 -O3
# FCOMP = ifort -r8 -i4 -O3
# FCOMP = xlf90 -qrealsize=8 -qintsize=4 -O3 -qstrict

all : rosseland

rosseland : RosselandOpacities.o
	$(FCOMP) -o $@ RosselandOpacities.o

.SUFFIXES: .F90

.F90.o:
	$(FCOMP) -c $*.F90

clean :
	rm -f *.o *.mod *~ rosseland

.PHONY: all
