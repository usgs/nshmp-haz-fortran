#!/bin/bash

time1=`date`
LOG=../logs
BIN=../bin
CONF=../conf

FAULT=$BIN/hazFXnga7c
GRID=$BIN/hazgridXnga5

CEUS_CONF=$CONF/CEUS
CEUS_FAULTS=1
CEUS_BACKGRND=0

# CEUS Faults
if [ $CEUS_FAULTS == 1 ]; then
	D_CEUS_F=$CEUS_CONF/faults
	$FAULT $D_CEUS_F/NMSZnocl.1000yr.5branch.in		> $LOG/nm_nocl1000y.log
	$FAULT $D_CEUS_F/NMSZnocl.500yr.5branch.in		> $LOG/nm_nocl500y.log
	$FAULT $D_CEUS_F/newmad.1000.cluster.in			> $LOG/newmad1000cl.log
	$FAULT $D_CEUS_F/newmad.1500.cluster.in			> $LOG/newmad1500cl.log
	$FAULT $D_CEUS_F/newmad.500.cluster.in			> $LOG/newmad500cl.log
	$FAULT $D_CEUS_F/newmad.750.cluster.in			> $LOG/newmad750cl.log
	$FAULT $D_CEUS_F/CEUScm.in 				> $LOG/ceuscm.log
fi

# CEUS Background 
if [ $CEUS_BACKGRND == 1 ]; then
	D_CEUS_BG=$CEUS_CONF/gridded
	$GRID $D_CEUS_BG/CEUS.2007all8.AB.in			> $LOG/ceusAB.log
	$GRID $D_CEUS_BG/CEUS.2007all8.J.in			> $LOG/ceusJ.log
	$GRID $D_CEUS_BG/CEUSchar.68.in				> $LOG/ceusChar68.log
	$GRID $D_CEUS_BG/CEUSchar.71.in			 	> $LOG/ceusChar71.log
	$GRID $D_CEUS_BG/CEUSchar.73.in				> $LOG/ceusChar73.log
	$GRID $D_CEUS_BG/CEUSchar.75.in				> $LOG/ceusChar75.log
fi
