# ============================================================================
# Name        : Makefile
# Author      : USGS
# Version     :
# Copyright   : 
# Description : Makefile for NSHMP14 codes
# ============================================================================

.PHONY: all clean

m_bit=

F_COMPILER = gfortran
F_COMPILER2 = ifort
FFLAGS1 = -O2 -Warray-bounds -ffixed-line-length-none -ffpe-trap= -fbounds-check $(m_bit) 
FFLAGS1_I = -132 
FFLAGS2 = $(FFLAGS1) $(m_bit)

C_COMPILER = gcc
CFLAGS = -O $(m_bit)

OUT = bin
SRC = src
UTIL = $(SRC)/util

## NOTE: Pick one or the other list of codes to compile (comment out the other!)
## compile codes on mac (no parallel version of hazFXnga131)
all: CreateBinDir hazallXL.v2 hazallXL.v4 hazFXnga13l hazgridXnga13l hazSUBX hazpoint hazinterpnga avg_dist fltrate.v2 get_avalue gethead.nga getmeanrjf getmeanrjf.v2 gutenberg assim.2013
## compile codes on CLUSTER (with parallel version of hazFXnga131)
#all: CreateBinDir hazallXL.v2 hazallXL.v4 hazFXnga13l hazFXnga13lp hazgridXnga13l hazSUBX hazpoint hazinterpnga avg_dist fltrate.v2 get_avalue gethead.nga getmeanrjf getmeanrjf.v2 gutenberg assim.2013

CreateBinDir:
	mkdir -p $(OUT)

##	dependencies
iosubs:
	$(C_COMPILER)  $(CFLAGS) -c -o $(SRC)/iosubs.o $(SRC)/iosubs.c
iosubs128:
	$(C_COMPILER)  $(CFLAGS) -c -o $(SRC)/iosubs_128.o $(SRC)/iosubs_128.c


##	hazard curve generation
hazallXL.v2: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/hazallXL.v2 $(SRC)/hazallXL.v2.f $(SRC)/iosubs.o
hazallXL.v4: iosubs
	$(F_COMPILER) $(FFLAGS2) -o $(OUT)/hazallXL.v4 $(SRC)/hazallXL.v4.f $(SRC)/iosubs.o 
	
##	fault hazard calculations (comment out hazFXnga131p if not compiling on CLUSTER!)
hazFXnga13l: iosubs
	$(F_COMPILER) $(FFLAGS1) -finit-local-zero -o $(OUT)/hazFXnga13l $(SRC)/hazFXnga13l.f $(SRC)/iosubs.o
#hazFXnga13lp: iosubs_128
	#$(F_COMPILER2) $(FFLAGS1_I) -o $(OUT)/hazFXnga13lp $(SRC)/hazFXnga13lp.f -coarray $(SRC)/iosubs_128.o -w

##	gridded hazard calculations
hazgridXnga13l: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/hazgridXnga13l $(SRC)/hazgridXnga13l.f $(SRC)/iosubs.o

##	ceus hazard calculations
hazSUBX: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/hazSUBX $(SRC)/hazSUBX.f $(SRC)/iosubs.o

## 	single-site hazard curve generation
hazpoint: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/hazpoint $(SRC)/hazpoint.f $(SRC)/iosubs.o
	
##	resampling	
hazinterpnga: iosubs
	$(F_COMPILER) $(FFLAGS2) -o $(OUT)/hazinterpnga $(SRC)/hazinterpnga.f $(SRC)/iosubs.o

## deagg codes
deaggFLTH: 
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/deaggFLTH $(SRC)/deaggFLTH.f 

deaggGRID: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/deaggGRID $(SRC)/deaggGRID.f $(SRC)/iosubs.o

deaggSUBD: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/deaggSUBD $(SRC)/deaggSUBD.f 

combine_cms: 
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/combine_cms $(SRC)/combine_cms.f 

sum_haz: 
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/sum_haz $(SRC)/sum_haz.f 

checksum09: 
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/checksum09 $(SRC)/checksum09.f 

##	utility
avg_dist: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/avg_dist $(UTIL)/avg_dist.f $(SRC)/iosubs.o
bin2xyzX.v2: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/bin2xyzX.v2 $(UTIL)/bin2xyzX.v2.f $(SRC)/iosubs.o
fltrate.v2: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/fltrate.v2 $(UTIL)/fltrate.v2.f $(SRC)/iosubs.o
fltrate.2013: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/fltrate.2013 $(UTIL)/fltrate.2013.f $(SRC)/iosubs.o
get_akprob: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/get_akprob $(UTIL)/get_akprob.f $(SRC)/iosubs.o
get_avalue: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/get_avalue $(UTIL)/get_avalue.f $(SRC)/iosubs.o
gethead.nga: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/gethead.nga $(UTIL)/gethead.nga.f $(SRC)/iosubs.o
gethead_noLong.nga: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/gethead_noLong.nga $(UTIL)/gethead.nga.f $(SRC)/iosubs_noLong.o
getmeanrjf: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/getmeanrjf $(UTIL)/getmeanrjf.f $(SRC)/iosubs.o
getmeanrjf.v2: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/getmeanrjf.v2 $(UTIL)/getmeanrjf.v2.f $(SRC)/iosubs.o
	scripts/meanrjb.sh
gutenberg: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/gutenberg $(UTIL)/gutenberg.f $(SRC)/iosubs.o
assim.2013: 
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/assim.2013 $(UTIL)/assim.2013.f
swapf: 
	$(C_COMPILER) $(CLAGS) -o $(OUT)/swapf $(UTIL)/swapf.c

clean:
	rm -f $(SRC)/*.o

cleanall:
	rm -f $(OUT)/* $(SRC)/*.o
