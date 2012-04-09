#!/bin/bash

time1=`date`
LOG=../logs
BIN=../bin
CONF=../conf

FAULT=$BIN/hazFXnga7c
GRID=$BIN/hazgridXnga5

CA_CONF=$CONF/CA
CA_A_FAULTS=1
CA_B_FAULTS=1
CA_BACKGRND=1
CA_SHEARZNS=1

# A Faults
if [ $CA_A_FAULTS == 1 ]; then
	D_CA_F=$CA_CONF/faults
	D_CA_F=$CA_CONF/faults
	$FAULT "$D_CA_F/aFault_aPriori_D2.1.in"			> $LOG/aF1.log 
#	$FAULT "$D_CA_F/aFault_MoBal_EllB.in"			> $LOG/aF2.log
#	$FAULT "$D_CA_F/aFault_MoBal.HB.in"			> $LOG/aF3.log 
#	$FAULT "$D_CA_F/aFault_unsegEll.in"			> $LOG/aF4.log 
#	$FAULT "$D_CA_F/aFault_unseg_HB.in"			> $LOG/aF5.log 
	$FAULT "$D_CA_F/aFault_unseg.in"			> $LOG/aF2.log 
	$FAULT "$D_CA_F/aFault_MoBal.in"			> $LOG/aF3.log 
#	$FAULT "$D_CA_F/creepflt.in"				> $LOG/creepflt.log 
fi


# B Faults
if [ $CA_B_FAULTS == 1 ]; then
	D_CA_F=$CA_CONF/faults
	$FAULT "$D_CA_F/bFault.ch.in"				> $LOG/bFaultch.log 
	$FAULT "$D_CA_F/bFault.gr.in"				> $LOG/bFaultgr.log 
fi


# Background Gridded Seismicity
if [ $CA_BACKGRND == 1 ]; then
	D_CA_BG=$CA_CONF/gridded
	$GRID "$D_CA_BG/CAdeep.in"				> $LOG/cadeep.log 
	$GRID "$D_CA_BG/CAmap.21.ch.in"				> $LOG/camap21ch.log 
	$GRID "$D_CA_BG/CAmap.24.ch.in"				> $LOG/camap24ch.log
	$GRID "$D_CA_BG/CAmap.21.gr.in"				> $LOG/camap21gr.log 
	$GRID "$D_CA_BG/CAmap.24.gr.in"				> $LOG/camap24gr.log
	$GRID "$D_CA_BG/brawmap.in"				> $LOG/brawmap.log 
	$GRID "$D_CA_BG/creepmap.in"				> $LOG/creepmmap.log 
fi


# Shear Zone Gridded Seismicity
if [ $CA_SHEARZNS == 1 ]; then
	D_CA_SZ=$CA_CONF/gridded
	$GRID "$D_CA_SZ/sangorg.in" 				> $LOG/sangorg.log 
	$GRID "$D_CA_SZ/mendo.in" 				> $LOG/mendo.log 
	$GRID "$D_CA_SZ/mojave.in"				> $LOG/mojave.log 
	$GRID "$D_CA_SZ/impext.ch.in"				> $LOG/impextch.log 
	$GRID "$D_CA_SZ/impext.gr.in"				> $LOG/impextgr.log 
	$GRID "$D_CA_SZ/shear1.in"				> $LOG/shear1.log 
	$GRID "$D_CA_SZ/shear2.in"				> $LOG/shear2.log 
	$GRID "$D_CA_SZ/shear3.in"				> $LOG/shear3.log 
	$GRID "$D_CA_SZ/shear4.in"				> $LOG/shear4.log 
fi

