#!/bin/bash

# Non-CA WUS fault hazard for 2014

# set local variables
time1=`date`
LOG=../logs
BIN=../bin
CONF=../conf
FAULT=$BIN/hazFXnga13l
GRID=$BIN/hazgridXnga13l
WUS_CONF=$CONF/WUS

WUS_FAULTS=1
WUS_GRIDDED=1

# Western US Faults
if (( $WUS_FAULTS == 1 )); then
	D_WUS_F=$WUS_CONF/faults
	
	# WUS faults - Bird geodetic
        $FAULT $D_WUS_F/2014WUSbird.65.in > $LOG/log_2014WUSbird.65.txt
        $FAULT $D_WUS_F/2014WUSbird.char.in > $LOG/log_2014WUSbird.char.txt
        $FAULT $D_WUS_F/2014WUSbird.gr.in > $LOG/log_2014WUSbird.gr.txt

	# WUS faults - Zeng geodetic
        $FAULT $D_WUS_F/2014WUSzeng.65.in > $LOG/log_2014WUSzeng.65.txt
        $FAULT $D_WUS_F/2014WUSzeng.char.in > $LOG/log_2014WUSzeng.char.txt
        $FAULT $D_WUS_F/2014WUSzeng.gr.in > $LOG/log_2014WUSzeng.gr.txt

	# WUS faults - geologic model
        $FAULT $D_WUS_F/2014WUSgeo.65.in > $LOG/log_2014WUSgeo.65.txt
        $FAULT $D_WUS_F/2014WUSgeo.char.in > $LOG/log_2014WUSgeo.char.txt
        $FAULT $D_WUS_F/2014WUSgeo.gr.in > $LOG/log_2014WUSgeo.gr.txt
        $FAULT $D_WUS_F/Seattle.char.in > $LOG/log_Seattle.char.txt
        $FAULT $D_WUS_F/Seattle.gr.in > $LOG/log_Seattle.gr.txt
        $FAULT $D_WUS_F/wasatch_slc.cluster.in > $LOG/log_wasatch_slc.cluster.txt
        $FAULT $D_WUS_F/wasatch_slc.noclu.in > $LOG/log_wasatch_slc.noclu.txt

fi


# Background Gridded Seismicity (7)
if (( $WUS_GRIDDED == 1 )); then
	D_WUS_G=$WUS_CONF/gridded
	
	# WUS-extensional
        $GRID $D_WUS_G/WUSmap_2014_adSm_r.ch.in > $LOG/log_WUSmap_2014_adSm_r.ch.txt
        $GRID $D_WUS_G/WUSmap_2014_adSm_r.gr.in > $LOG/log_WUSmap_2014_adSm_r.gr.txt
        $GRID $D_WUS_G/WUSmap_2014_adSm_r_M8.in > $LOG/log_WUSmap_2014_adSm_r_M8.txt
        $GRID $D_WUS_G/WUSmap_2014_adSm_ss.ch.in > $LOG/log_WUSmap_2014_adSm_ss.ch.txt
        $GRID $D_WUS_G/WUSmap_2014_adSm_ss.gr.in > $LOG/log_WUSmap_2014_adSm_ss.gr.txt
        $GRID $D_WUS_G/WUSmap_2014_adSm_ss_M8.in > $LOG/log_WUSmap_2014_adSm_ss_M8.txt
        $GRID $D_WUS_G/WUSmap_2014_fixSm_r.ch.in > $LOG/log_WUSmap_2014_fixSm_r.ch.txt
        $GRID $D_WUS_G/WUSmap_2014_fixSm_r.gr.in > $LOG/log_WUSmap_2014_fixSm_r.gr.txt
        $GRID $D_WUS_G/WUSmap_2014_fixSm_r_M8.in > $LOG/log_WUSmap_2014_fixSm_r_M8.txt
        $GRID $D_WUS_G/WUSmap_2014_fixSm_ss.ch.in > $LOG/log_WUSmap_2014_fixSm_ss.ch.txt
        $GRID $D_WUS_G/WUSmap_2014_fixSm_ss.gr.in > $LOG/log_WUSmap_2014_fixSm_ss.gr.txt
        $GRID $D_WUS_G/WUSmap_2014_fixSm_ss_M8.in > $LOG/log_WUSmap_2014_fixSm_ss_M8.txt


        # WUS-compressional
        $GRID $D_WUS_G/WUSmap_2014_adSm_r.ch.in > $LOG/log_WUSmap_2014_adSm_r.ch.txt
        $GRID $D_WUS_G/WUSmap_2014_adSm_r.gr.in > $LOG/log_WUSmap_2014_adSm_r.gr.txt
        $GRID $D_WUS_G/WUSmap_2014_adSm_r_M8.in > $LOG/log_WUSmap_2014_adSm_r_M8.txt
        $GRID $D_WUS_G/WUSmap_2014_adSm_ss.ch.in > $LOG/log_WUSmap_2014_adSm_ss.ch.txt
        $GRID $D_WUS_G/WUSmap_2014_adSm_ss.gr.in > $LOG/log_WUSmap_2014_adSm_ss.gr.txt
        $GRID $D_WUS_G/WUSmap_2014_adSm_ss_M8.in > $LOG/log_WUSmap_2014_adSm_ss_M8.txt
        $GRID $D_WUS_G/WUSmap_2014_fixSm_r.ch.in > $LOG/log_WUSmap_2014_fixSm_r.ch.txt
        $GRID $D_WUS_G/WUSmap_2014_fixSm_r.gr.in > $LOG/log_WUSmap_2014_fixSm_r.gr.txt
        $GRID $D_WUS_G/WUSmap_2014_fixSm_r_M8.in > $LOG/log_WUSmap_2014_fixSm_r_M8.txt
        $GRID $D_WUS_G/WUSmap_2014_fixSm_ss.ch.in > $LOG/log_WUSmap_2014_fixSm_ss.ch.txt
        $GRID $D_WUS_G/WUSmap_2014_fixSm_ss.gr.in > $LOG/log_WUSmap_2014_fixSm_ss.gr.txt
        $GRID $D_WUS_G/WUSmap_2014_fixSm_ss_M8.in > $LOG/log_WUSmap_2014_fixSm_ss_M8.txt
        # WUS-compressional, complement Puget region
        $GRID $D_WUS_G/noPuget_2014_adSm_r.ch.in > $LOG/log_noPuget_2014_adSm_r.ch.txt
        $GRID $D_WUS_G/noPuget_2014_adSm_r.gr.in > $LOG/log_noPuget_2014_adSm_r.gr.txt
        $GRID $D_WUS_G/noPuget_2014_adSm_r_M8.in > $LOG/log_noPuget_2014_adSm_r_M8.txt
        $GRID $D_WUS_G/noPuget_2014_adSm_ss.ch.in > $LOG/log_noPuget_2014_adSm_ss.ch.txt
        $GRID $D_WUS_G/noPuget_2014_adSm_ss.gr.in > $LOG/log_noPuget_2014_adSm_ss.gr.txt
        $GRID $D_WUS_G/noPuget_2014_adSm_ss_M8.in > $LOG/log_noPuget_2014_adSm_ss_M8.txt
        $GRID $D_WUS_G/noPuget_2014_fixSm_r.ch.in > $LOG/log_noPuget_2014_fixSm_r.ch.txt
        $GRID $D_WUS_G/noPuget_2014_fixSm_r.gr.in > $LOG/log_noPuget_2014_fixSm_r.gr.txt
        $GRID $D_WUS_G/noPuget_2014_fixSm_r_M8.in > $LOG/log_noPuget_2014_fixSm_r_M8.txt
        $GRID $D_WUS_G/noPuget_2014_fixSm_ss.ch.in > $LOG/log_noPuget_2014_fixSm_ss.ch.txt
        $GRID $D_WUS_G/noPuget_2014_fixSm_ss.gr.in > $LOG/log_noPuget_2014_fixSm_ss.gr.txt
        $GRID $D_WUS_G/noPuget_2014_fixSm_ss_M8.in > $LOG/log_noPuget_2014_fixSm_ss_M8.txt
        $GRID $D_WUS_G/puget_2014_r.ch.in > $LOG/log_puget_2014_r.ch.txt
        $GRID $D_WUS_G/puget_2014_r.gr.in > $LOG/log_puget_2014_r.gr.txt
        $GRID $D_WUS_G/puget_2014_r_M8.in > $LOG/log_puget_2014_r_M8.txt
        $GRID $D_WUS_G/puget_2014_ss.ch.in > $LOG/log_puget_2014_ss.ch.txt
        $GRID $D_WUS_G/puget_2014_ss.gr.in > $LOG/log_puget_2014_ss.gr.txt
        $GRID $D_WUS_G/puget_2014_ss_M8.in > $LOG/log_puget_2014_ss_M8.txt


        # gridded fault hazard from Bird geodetic model
        $GRID $D_WUS_G/WUS_zones_PB.in > $LOG/log_WUS_zones_PB.txt

        # WUS shear zones
        $GRID $D_WUS_G/shear2_2014.in > $LOG/log_shear2_2014.txt
        $GRID $D_WUS_G/shear3_2014.in > $LOG/log_shear3_2014.txt
        $GRID $D_WUS_G/shear4_2014.in > $LOG/log_shear4_2014.txt

        # deep seismicity - CA, coastal OR, puget sound
        $GRID $D_WUS_G/CAdeep.2014.in > $LOG/log_CAdeep.2014.txt
        $GRID $D_WUS_G/CAdeepMmax75.2014.in > $LOG/log_CAdeepMmax75.2014.txt
        $GRID $D_WUS_G/CAdeepMmax8.2014.in > $LOG/log_CAdeepMmax8.2014.txt
        $GRID $D_WUS_G/coastalOR_deep.in > $LOG/log_coastalOR_deep.txt
        $GRID $D_WUS_G/coastalOR_deep_Mmax75.in > $LOG/log_coastalOR_deep_Mmax75.txt
        $GRID $D_WUS_G/coastalOR_deep_Mmax8.in > $LOG/log_coastalOR_deep_Mmax8.txt
        $GRID $D_WUS_G/pacnwdeep.2014.in > $LOG/log_pacnwdeep.2014.txt
        $GRID $D_WUS_G/pacnwdeep_Mmax75.2014.in > $LOG/log_pacnwdeep_Mmax75.2014.txt
        $GRID $D_WUS_G/pacnwdeep_Mmax8.2014.in > $LOG/log_pacnwdeep_Mmax8.2014.txt

fi

echo "$0 script completed"

