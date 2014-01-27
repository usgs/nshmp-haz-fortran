# ============================================================================
# Name        : Makefile
# Author      : USGS
# Version     :
# Copyright   : 
# Description : Makefile for NSHMP08 codes
# ============================================================================

.PHONY: all clean

#m_bit=-m32 # flag does not work on Linux cluster
#m_bit=-m32 # flag does not work on Linux cluster
m_bit=

F_COMPILER = gfortran
F_COMPILER2 = ifort
FFLAGS1 = -O2 -Warray-bounds -ffixed-line-length-none -ffpe-trap= -fbounds-check $(m_bit) 
FFLAGS1_I = -132 
#FFLAGS2 = $(FFLAGS1) -fcray-pointer $(m_bit)
FFLAGS2 = $(FFLAGS1) $(m_bit)

C_COMPILER = gcc
CFLAGS = -O $(m_bit)

OUT = bin
SRC = src
UTIL = $(SRC)/util


#all: CreateBinDir hazallXL.v2 hazallXL.v4 hazFXnga7c hazFXnga13l hazgridXnga5 hazgridXnga13l hazSUBXnga hazSUBXngatest hazpoint hazinterpnga avg_dist fltrate.v2 get_avalue gethead.nga getmeanrjf getmeanrjf.v2 gutenberg
#all: CreateBinDir hazallXL.v2 hazallXL.v4 hazallXL.v5 hazFXnga13l hazgridXnga13l hazSUBXnga hazSUBXngatest hazpoint hazinterpnga avg_dist fltrate.v2 get_avalue gethead.nga getmeanrjf getmeanrjf.v2 gutenberg
#all: CreateBinDir hazallXL.v2 hazallXL.v4 hazFXnga13l hazgridXnga13l hazSUBXnga hazSUBXngatest hazpoint hazinterpnga avg_dist fltrate.v2 get_avalue gethead.nga getmeanrjf getmeanrjf.v2 gutenberg
all: CreateBinDir hazallXL.v2 hazallXL.v4 hazFXnga13l hazgridXnga13l hazSUBX hazpoint hazinterpnga avg_dist fltrate.v2 get_avalue gethead.nga getmeanrjf getmeanrjf.v2 gutenberg assim.2013

CreateBinDir:
	mkdir -p $(OUT)

#	dependencies
iosubs:
	$(C_COMPILER)  $(CFLAGS) -c -o $(SRC)/iosubs.o $(SRC)/iosubs.c
iosubs128:
	$(C_COMPILER)  $(CFLAGS) -c -o $(SRC)/iosubs_128.o $(SRC)/iosubs_128.c
iosubs_noLong:
	$(C_COMPILER)  $(CFLAGS) -c -o $(SRC)/iosubs_noLong.o $(SRC)/iosubs_noLong.c

#	hazard curve generation
hazallXL.v2: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/hazallXL.v2 $(SRC)/hazallXL.v2.f $(SRC)/iosubs.o
hazallXL.v4: iosubs
	$(F_COMPILER) $(FFLAGS2) -o $(OUT)/hazallXL.v4 $(SRC)/hazallXL.v4.f $(SRC)/iosubs.o 
hazallXL.v5: iosubs
	$(F_COMPILER) $(FFLAGS2) -o $(OUT)/hazallXL.v5 $(SRC)/hazallXL.v5.f $(SRC)/iosubs.o 
#hazallXL.v5: iosubs128
#	$(F_COMPILER) $(FFLAGS2) -o $(OUT)/hazallXL.v5 $(SRC)/hazallXL.v5.f $(SRC)/iosubs_128.o 

hazFXnga7c: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/hazFXnga7c $(SRC)/hazFXnga7c.f $(SRC)/iosubs.o
hazFXnga13l: iosubs
	$(F_COMPILER) $(FFLAGS1) -finit-local-zero -o $(OUT)/hazFXnga13l $(SRC)/hazFXnga13l.f $(SRC)/iosubs.o
hazFXnga13p: iosubs128
	ifort -o $(OUT)/hazFXnga13p $(SRC)/hazFXnga13p.f -coarray $(SRC)/iosubs_128.o -w -132
hazFXnga7.temp: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/hazFXnga7.temp $(SRC)/hazFXnga7.temp.f $(SRC)/iosubs.o
hazFXnga12: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/hazFXnga12 $(SRC)/hazFXnga12.f $(SRC)/iosubs.o

hazgridXnga5: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/hazgridXnga5 $(SRC)/hazgridXnga5.f $(SRC)/iosubs.o
hazgridXnga13l: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/hazgridXnga13l $(SRC)/hazgridXnga13l.f $(SRC)/iosubs.o
hazgridXnga13lidOld: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/hazgridXnga13lidOld $(SRC)/hazgridXnga13lidOld.f $(SRC)/iosubs.o
hazgridXnga13l_test: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/hazgridXnga13l_test $(SRC)/hazgridXnga13l_test.f $(SRC)/iosubs.o
hazgridXnga13l_deep: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/hazgridXnga13l_deep $(SRC)/hazgridXnga13l_deep.f $(SRC)/iosubs.o

hazSUBXnga: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/hazSUBXnga $(SRC)/hazSUBXnga.f $(SRC)/iosubs.o
hazSUBXngatest: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/hazSUBXngatest $(SRC)/hazSUBXngatest.f $(SRC)/iosubs.o
hazSUBX: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/hazSUBX $(SRC)/hazSUBX.f $(SRC)/iosubs.o

hazpoint: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/hazpoint $(SRC)/hazpoint.f $(SRC)/iosubs.o
hazinterpnga: iosubs
#	$(F_COMPILER2) $(FFLAGS1_I) -o $(OUT)/hazinterpnga $(SRC)/hazinterpnga.f $(SRC)/iosubs_128.o -w
	$(F_COMPILER) $(FFLAGS2) -o $(OUT)/hazinterpnga $(SRC)/hazinterpnga.f $(SRC)/iosubs.o

# deagg codes
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

#	utility
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
gutenberg: iosubs
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/gutenberg $(UTIL)/gutenberg.f $(SRC)/iosubs.o
assim.2013: 
	$(F_COMPILER) $(FFLAGS1) -o $(OUT)/assim.2013 $(UTIL)/assim.2013.f
swapf: 
	$(C_COMPILER) $(CLAGS) -o $(OUT)/swapf $(UTIL)/swapf.c

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
