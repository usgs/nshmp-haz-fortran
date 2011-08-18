#!/bin/bash

time1=`date`

LOG=../logs
BIN=../bin
CONF=../conf

SUB1=$BIN/hazSUBXnga
SUB2=$BIN/hazSUBXngatest

CASC_CONF=$CONF/CASC
CASC_SUB=1

# All Subduction
if [ $CASC_SUB == 1 ]; then

	$SUB2 "$CASC_CONF/cascadia.bot.9pm.in"		> $LOG/cascBot9.log 
	$SUB2 "$CASC_CONF/cascadia.mid.9pm.in"		> $LOG/cascMid9.log 
	$SUB2 "$CASC_CONF/cascadia.top.9pm.in"		> $LOG/cascTop9.log 
	$SUB2 "$CASC_CONF/cascadia.older2.9pm.in"	> $LOG/cascOld9.log 
	
	$SUB1 "$CASC_CONF/cascadia.bot.8082.in"		> $LOG/cascBot8082.log 
	$SUB1 "$CASC_CONF/cascadia.mid.8082.in"		> $LOG/cascMid8082.log 
	$SUB1 "$CASC_CONF/cascadia.top.8082.in"		> $LOG/cascTop8082.log 
	$SUB1 "$CASC_CONF/cascadia.older2.8082.in"	> $LOG/cascOld8082.log 
	$SUB1 "$CASC_CONF/cascadia.bot.8387.in"		> $LOG/cascBot8387.log 
	$SUB1 "$CASC_CONF/cascadia.mid.8387.in"		> $LOG/cascMid8387.log 
	$SUB1 "$CASC_CONF/cascadia.top.8387.in"		> $LOG/cascTop8387.log 
	$SUB1 "$CASC_CONF/cascadia.older2.8387.in"	> $LOG/cascOld8387.log 
fi

