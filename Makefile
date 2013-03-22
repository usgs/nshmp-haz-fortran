# ============================================================================
# Name        : Makefile
# Author      : USGS
# Version     :
# Copyright   : 
# Description : Makefile for NSHMP08 codes
# ============================================================================

.PHONY: all clean

#m_bit=-m32
m_bit=-m32

F_COMPILER = gfortran
FFLAGS1 = -O2 -Warray-bounds -ffixed-line-length-none -ffpe-trap= $(m_bit)
FFLAGS2 = $(FFLAGS1) -fcray-pointer $(m_bit)

C_COMPILER = gcc
CFLAGS = -O $(m_bit)

OUT = bin
SRC = src
UTIL = $(SRC)/util


all:
	mkdir -p $(OUT)

#	dependencies
	$(C_COMPILER)  $(CFLAGS) -c -o $(SRC)/iosubs.o $(SRC)/iosubs.c
	$(C_COMPILER)  $(CFLAGS) -c -o $(SRC)/iosubs_128.o $(SRC)/iosubs_128.c

#	hazard curve generation
	$(F_COMPILER) $(FFLAGS2) -o $(OUT)/hazallXL.v4 $(SRC)/hazallXL.v4.f $(SRC)/iosubs.o 
	$(F_COMPILER) $(FFLAGS2) -o $(OUT)/hazallXL.v5 $(SRC)/hazallXL.v5.f $(SRC)/iosubs_128.o 

	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/hazFXnga7c $(SRC)/hazFXnga7c.f $(SRC)/iosubs.o
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/hazFXnga13l $(SRC)/hazFXnga13l.f $(SRC)/iosubs.o
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/hazFXnga7.temp $(SRC)/hazFXnga7.temp.f $(SRC)/iosubs.o
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/hazFXnga12 $(SRC)/hazFXnga12.f $(SRC)/iosubs.o
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/hazgridXnga5 $(SRC)/hazgridXnga5.f $(SRC)/iosubs.o
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/hazgridXnga13l $(SRC)/hazgridXnga13l.f $(SRC)/iosubs.o
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
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/getmeanrjf.v2 $(UTIL)/getmeanrjf.v2.f $(SRC)/iosubs.o
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
