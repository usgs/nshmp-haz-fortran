#!/bin/bash

# Non-CA WUS fault hazard for 2008   plus Goose Lake in CA
# incorporates file consolidation of dip-uncertainty branches
# Basin and Range states AZ, MT, ID, CO, TX, WY, and NM
# Oregon is now divided into two models, coastal compressional and BR extensional
# orwa_c.in contains sources of orwa_c.char and seattle fault
# new oct 4 basin and range portion of Oregon. These have different wts in the logic tree
# compressional-src wts are 50% char 50% GR, extensional are 67%char 33% GR.

# orwa_c.in includes:
#	1) Both char and gr models for most compressional WUS faults
#		NOTE: there a few offshore faults with very low wight (e.g. Alvin
#		Canyon, Wecoma, and Daisy Bank @ 0.05 wieght) 
#	2) M < 6.5 events (was orwa.65) with double weight because they're
#		modeled as characteristic only.
#	3) Seattle Fault Zone (was SeattleFZ.in)
#		NOTE: By virtue of this combination, orwa_c.ch.in includes 
#		the GR components of the Seattle Fault. This is ok for now 
#		as the GR and CH parts of compressional faults in the PNW
#		each get 50% weight.

# Input files for normal faults contain sources with 40, 50, and 60 degree
# dips with 0.2, 0.6, and 0.2 weights respectively. These files used to have
# '3dip' as part of their name
# the floating M7.4 along wasatch, slip rate is 1.2 mm/yr for this model MP Feb 9 2007

# set local variables
time1=`date`
LOG=../logs
BIN=../bin
CONF=../conf
FAULT=$BIN/hazFXnga7c
GRID=$BIN/hazgridXnga5
WUS_CONF=$CONF/WUS

WUS_FAULTS=1
WUS_GRIDDED=1

# Western US Faults
if (( $WUS_FAULTS == 1 )); then
	D_WUS_F=$WUS_CONF/faults
	
	# Basin and Range
	$FAULT "$D_WUS_F/brange.3dip.ch.in"		> $LOG/brangech.log 
	$FAULT "$D_WUS_F/brange.3dip.gr.in" 	> $LOG/brangegr.log
	$FAULT "$D_WUS_F/brange.3dip.65.in" 	> $LOG/brange65.log 

	# PNW: Oregon and Washington
	$FAULT "$D_WUS_F/orwa_c.in" 			> $LOG/orwac.log 
	$FAULT "$D_WUS_F/orwa_n.3dip.ch.in" 	> $LOG/orwanch.log 
	$FAULT "$D_WUS_F/orwa_n.3dip.gr.in" 	> $LOG/orwangr.log

	# Nevada
	$FAULT "$D_WUS_F/nv.3dip.ch.in" 		> $LOG/nvch.log 
	$FAULT "$D_WUS_F/nv.3dip.gr.in" 		> $LOG/nvgr.log
	$FAULT "$D_WUS_F/nvut.3dip.65.in" 		> $LOG/nvut65.log 

	# Utah except Wasatch
	$FAULT "$D_WUS_F/ut.3dip.ch.in" 		> $LOG/utch.log 
	$FAULT "$D_WUS_F/ut.3dip.gr.in" 		> $LOG/utgr.log

	# Wasatch
	$FAULT "$D_WUS_F/wasatch.3dip.ch.in" 	> $LOG/wasatchch.log 
	$FAULT "$D_WUS_F/wasatch.3dip.gr.in" 	> $LOG/wasatchgr.log
	$FAULT "$D_WUS_F/wasatch.3dip.74.in"	> $LOG/wasatch74.log 

fi


# Background Gridded Seismicity (7)
if (( $WUS_GRIDDED == 1 )); then
	D_WUS_G=$WUS_CONF/gridded
	
	$GRID "$D_WUS_G/WUSmap.ch.in"			> $LOG/wusmapch.log 
	$GRID "$D_WUS_G/WUSmap.gr.in"			> $LOG/wusmapgr.log 
	$GRID "$D_WUS_G/EXTmap.ch.in"			> $LOG/extmapch.log 
	$GRID "$D_WUS_G/EXTmap.gr.in"			> $LOG/extmapgr.log 
	# Pacific Northwest
	$GRID "$D_WUS_G/nopuget.ch.in"			> $LOG/nopugetch.log 
	$GRID "$D_WUS_G/nopuget.gr.in"			> $LOG/nopugetgr.log 
	$GRID "$D_WUS_G/puget.ch.in"			> $LOG/pugetch.log 
	$GRID "$D_WUS_G/puget.gr.in"			> $LOG/pugetgr.log 
	$GRID "$D_WUS_G/pnwdeep.in"			> $LOG/pnwdeep.log 
	$GRID "$D_WUS_G/portdeep.in"			> $LOG/portdeep.log 

fi

echo "$0 script completed"

