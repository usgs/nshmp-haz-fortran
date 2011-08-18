# ============================================================================
# Name        : Makefile
# Author      : USGS
# Version     :
# Copyright   : 
# Description : Makefile for NSHMP08 codes
# ============================================================================

.PHONY: all clean

FORTRAN_COMPILER = gfortran
FFLAGS = -O2 -Warray-bounds -ffixed-line-length-none -ffpe-trap= 
#C_COMPILER = cc
C_COMPILER = gcc
CFLAGS = -O -arch i386 

OUT = bin
SRC = src

#
#all: src/NSHMP08.f90
#	$(FORTRAN_COMPILER) -O2 -g \
#		-o bin/NSHMP08 \
#		src/NSHMP08.f90


all:
#	gfortran hazFXnga7c.f -ffixed-line-length-none -ffpe-trap= -o hazFXnga7c -finit-local-zero
#	gfortran -ffixed-line-length-none -ffpe-trap= -o hazFXnga7c -g hazFXnga7c.f iosubs.o
#	cc -O -arch i386 -c swapf.c
#	gfortran -ffixed-line-length-none -fbounds-check -ffpe-trap= -o getmeanrjf getmeanrjf.f
#	cc -O -arch i386 -o swapf swapf.c
#	-fbounds-check 

	$(C_COMPILER)  $(CFLAGS) -c -o $(SRC)/iosubs.o $(SRC)/iosubs.c
	$(FORTRAN_COMPILER) $(FFLAGS) -fcray-pointer -o $(OUT)/hazallXL.v4 $(SRC)/hazallXL.v4.f $(SRC)/iosubs.o 
	$(FORTRAN_COMPILER) $(FFLAGS) -o $(OUT)/hazFXnga7c $(SRC)/hazFXnga7c.f $(SRC)/iosubs.o
	$(FORTRAN_COMPILER) $(FFLAGS) -o $(OUT)/hazgridXnga5 $(SRC)/hazgridXnga5.f $(SRC)/iosubs.o
	$(FORTRAN_COMPILER) $(FFLAGS) -o $(OUT)/hazSUBXnga $(SRC)/hazSUBXnga.f $(SRC)/iosubs.o
	$(FORTRAN_COMPILER) $(FFLAGS) -o $(OUT)/hazSUBXngatest $(SRC)/hazSUBXngaTest.f $(SRC)/iosubs.o
	$(FORTRAN_COMPILER) $(FFLAGS) -o $(OUT)/hazallXLv2 $(SRC)/hazallXLv2.f $(SRC)/iosubs.o
	$(FORTRAN_COMPILER) $(FFLAGS) -o $(OUT)/gethead.nga $(SRC)/gethead.nga.f $(SRC)/iosubs.o
	$(FORTRAN_COMPILER) $(FFLAGS) -fcray-pointer -o $(OUT)/hazinterpnga $(SRC)/hazinterpnga.f $(SRC)/iosubs.o
	$(FORTRAN_COMPILER) $(FFLAGS) -o $(OUT)/getmeanrjf $(SRC)/getmeanrjf.f $(SRC)/iosubs.o

	
clean:
	rm -f $(SRC)/*.o

cleanall:
	rm -f $(OUT)/* $(SRC)/*.o
