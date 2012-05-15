# ============================================================================
# Name        : Makefile
# Author      : USGS
# Version     :
# Copyright   : 
# Description : Makefile for NSHMP08 codes
# ============================================================================

.PHONY: all clean

F_COMPILER = gfortran
FFLAGS1 = -O2 -Warray-bounds -ffixed-line-length-none -ffpe-trap=
FFLAGS2 = $(FFLAGS1) -fcray-pointer

C_COMPILER = gcc
CFLAGS = -O -arch i386 

OUT = bin
SRC = src
UTIL = $(SRC)/util


all:
	mkdir -p $(OUT)

#	dependencies
	$(C_COMPILER)  $(CFLAGS) -c -o $(SRC)/iosubs.o $(SRC)/iosubs.c

#	hazard curve generation
	$(F_COMPILER) $(FFLAGS2) -o $(OUT)/hazallXL.v4 $(SRC)/hazallXL.v4.f $(SRC)/iosubs.o 
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/hazFXnga7c $(SRC)/hazFXnga7c.f $(SRC)/iosubs.o
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/hazgridXnga5 $(SRC)/hazgridXnga5.f $(SRC)/iosubs.o
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/hazSUBXnga $(SRC)/hazSUBXnga.f $(SRC)/iosubs.o
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/hazSUBXngatest $(SRC)/hazSUBXngaTest.f $(SRC)/iosubs.o
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/hazallXLv2 $(SRC)/hazallXLv2.f $(SRC)/iosubs.o
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/hazpoint $(SRC)/hazpoint.f $(SRC)/iosubs.o
	$(F_COMPILER) $(FFLAGS2) -o $(OUT)/hazinterpnga $(SRC)/hazinterpnga.f $(SRC)/iosubs.o

#	utility
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/avg_dist $(UTIL)/avg_dist.f $(SRC)/iosubs.o
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/fltrate.v2 $(UTIL)/fltrate.v2.f $(SRC)/iosubs.o
#	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/get_akprob $(UTIL)/get_akprob.f $(SRC)/iosubs.o
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/get_avalue $(UTIL)/get_avalue.f $(SRC)/iosubs.o
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/gethead.nga $(UTIL)/gethead.nga.f $(SRC)/iosubs.o
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/getmeanrjf $(UTIL)/getmeanrjf.f $(SRC)/iosubs.o
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/gutenberg $(UTIL)/gutenberg.f $(SRC)/iosubs.o

#	other
#	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/agridtest2A $(SRC)/agridtest2A.f 
#
#	gfortran hazFXnga7c.f -ffixed-line-length-none -ffpe-trap= -o hazFXnga7c -finit-local-zero
#	gfortran -ffixed-line-length-none -ffpe-trap= -o hazFXnga7c -g hazFXnga7c.f iosubs.o
#	cc -O -arch i386 -c swapf.c
#	gfortran -ffixed-line-length-none -fbounds-check -ffpe-trap= -o getmeanrjf getmeanrjf.f
#	cc -O -arch i386 -o swapf swapf.c
#	-fbounds-check 


clean:
	rm -f $(SRC)/*.o

cleanall:
	rm -f $(OUT)/* $(SRC)/*.o
