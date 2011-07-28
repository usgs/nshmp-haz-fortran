CC=gcc -Wall -O3 -static
FF=gfortran -Wall -O3 -static -ffixed-line-length-none


EXECUTABLES=\
	bin2xyzX.v2 \
	checksum09 \
	checksum6 \
	checksum \
	combine.2008b \
	combine.2008 \
	combine.2009 \
	combinehaz.2008b \
	craton_ck \
	deaggFLTH.2008 \
	deaggFLTH.2009 \
	deaggFLTH.PGA08 \
	deaggGRID.2008 \
	deaggGRID.2009 \
	deaggGRID.PGA08 \
	deaggSUBD.2008 \
	deaggSUBD.2009 \
	deaggSUBD.PGA08 \
	decide6 \
	decides \
	inwus \
	rescale \
	rjbmean_ck \
	sum_haz_2008b \
	sum_haz_2008 \
	sum_haz_2009 \
	swapf \
	vancouver

OBJECTS=\
	datetime.o \
	iosubs.o \
	recipes.o \
	rvtdsubs.o \
	smsim_tdww2.o

#####################################
# Static file must exist "smsim.fi" #
#####################################

#-------------------------------------------------------------------------------
# Generic Targets
#-------------------------------------------------------------------------------
install: objects exec
	cp src/run_2009.csh run_2009
	cp src/install.sh install
	chmod +x install
	./install
	rm -f install
	mkdir LOG
	@echo 'Build completed successfully. See above for any warnings.'

objects: init $(OBJECTS)

exec: objects $(EXECUTABLES)

init:
	-mkdir bin

clean:
	-rm -rf LOG
	-rm -rf bin
	-rm -f run_2009
	cp src/uninstall.sh uninstall
	chmod +x uninstall
	./uninstall
	-rm -f uninstall

#-------------------------------------------------------------------------------
# Object Targets
#-------------------------------------------------------------------------------

datetime.o:
	$(FF) -c src/datetime.f -o bin/datetime.o 
	
iosubs.o:
	$(CC) -c src/iosubs.c -o bin/iosubs.o

recipes.o:
	$(FF) -c src/recipes.f -o bin/recipes.o

rvtdsubs.o:
	$(FF) -c src/rvtdsubs.f -o bin/rvtdsubs.o

smsim_tdww2.o:
	$(FF) -c src/smsim_tdww2.f -o bin/smsim_tdww2.o

#-------------------------------------------------------------------------------
# Executable Targets
#-------------------------------------------------------------------------------

bin2xyzX.v2: iosubs.o
	$(FF) bin/iosubs.o -o bin/$@ src/$@.f

checksum09:
	$(FF) -o bin/$@ src/$@.f

checksum6:
	$(FF) -o bin/$@ src/$@.f

checksum:
	$(FF) -o bin/$@ src/$@.f

combine.2008b: iosubs.o
	$(FF) bin/iosubs.o -o bin/$@ src/$@.f

combine.2008: iosubs.o
	$(FF) bin/iosubs.o -o bin/$@ src/$@.f

combine.2009:iosubs.o
	$(FF) bin/iosubs.o -o bin/$@ src/$@.f

combinehaz.2008b: iosubs.o
	$(FF) bin/iosubs.o -o bin/$@ src/$@.f

craton_ck:
	$(FF) -o bin/$@ src/$@.f

deaggFLTH.2008: iosubs.o
	$(FF) bin/iosubs.o -o bin/$@ src/$@.f

deaggFLTH.2009: iosubs.o
	$(FF) bin/iosubs.o -o bin/$@ src/$@.f

deaggFLTH.PGA08: iosubs.o
	$(FF) bin/iosubs.o -o bin/$@ src/$@.f

deaggGRID.2008: iosubs.o
	$(FF) bin/iosubs.o -o bin/$@ src/$@.f

deaggGRID.2009: iosubs.o
	$(FF) bin/iosubs.o -o bin/$@ src/$@.f

deaggGRID.PGA08: iosubs.o
	$(FF) bin/iosubs.o -o bin/$@ src/$@.f

deaggSUBD.2008: iosubs.o
	$(FF) bin/iosubs.o -o bin/$@ src/$@.f

deaggSUBD.2009: iosubs.o
	$(FF) bin/iosubs.o -o bin/$@ src/$@.f

deaggSUBD.PGA08: iosubs.o
	$(FF) bin/iosubs.o -o bin/$@ src/$@.f

decide6:
	$(FF) -o bin/$@ src/$@.f

decides:
	$(FF) -o bin/$@ src/$@.f

inwus:
	$(FF) -o bin/$@ src/$@.f

rescale:
	$(FF) -o bin/$@ src/$@.f

rjbmean_ck:
	$(FF) -o bin/$@ src/$@.f

smsim.fi:
	@echo 'Nothing to be done for smsim.fi'

sum_haz_2008b:
	$(FF) -o bin/$@ src/$@.f

sum_haz_2008: iosubs.o recipes.o smsim_tdww2.o
	$(FF) bin/datetime.o bin/iosubs.o bin/recipes.o bin/rvtdsubs.o bin/smsim_tdww2.o -o bin/$@ src/$@.f

sum_haz_2009: iosubs.o recipes.o smsim_tdww2.o
	$(FF) bin/datetime.o bin/iosubs.o bin/recipes.o bin/rvtdsubs.o bin/smsim_tdww2.o -o bin/$@ src/$@.f

swapf:
	$(CC) -o bin/$@ src/$@.c

vancouver:
	$(FF) -o bin/$@ src/$@.f
