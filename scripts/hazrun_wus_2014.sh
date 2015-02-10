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

WUS_FAULTS=0
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
        for inputf in 2014WUSgeo.65.in 2014WUSgeo.char.in 2014WUSgeo.gr.in Seattle.char.in Seattle.gr.in wasatch_slc.cluster.in wasatch_slc.noclu.in ; do
          logf=`echo $inputf | sed 's/in$/txt/' | awk '{print "../logs/log_"$1}'`
          echo "Running $inputf..."
          echo "logfile: $logf"
          $FAULT $D_WUS_F/$inputf > $logf
        done

fi


# Background Gridded Seismicity (7)
if (( $WUS_GRIDDED == 1 )); then
	D_WUS_G=$WUS_CONF/gridded
	
	# WUS-extensional
	for inputf in EXTmap_2014_adSm_n.ch.in EXTmap_2014_adSm_n.gr.in EXTmap_2014_adSm_n_M8.in EXTmap_2014_adSm_ss.ch.in EXTmap_2014_adSm_ss.gr.in EXTmap_2014_adSm_ss_M8.in EXTmap_2014_fixSm_n.ch.in EXTmap_2014_fixSm_n.gr.in EXTmap_2014_fixSm_n_M8.in EXTmap_2014_fixSm_ss.ch.in EXTmap_2014_fixSm_ss.gr.in EXTmap_2014_fixSm_ss_M8.in; do
          logf=`echo $inputf | sed 's/in$/txt/' | awk '{print "../logs/log_"$1}'`
          echo "Running $inputf..."
          echo "logfile: $logf"
          $GRID $D_WUS_G/$inputf > $logf
        done

        # WUS-compressional
	for inputf in WUSmap_2014_adSm_r.ch.in WUSmap_2014_adSm_r.gr.in WUSmap_2014_adSm_r_M8.in WUSmap_2014_adSm_ss.ch.in WUSmap_2014_adSm_ss.gr.in WUSmap_2014_adSm_ss_M8.in WUSmap_2014_fixSm_r.ch.in WUSmap_2014_fixSm_r.gr.in WUSmap_2014_fixSm_r_M8.in WUSmap_2014_fixSm_ss.ch.in WUSmap_2014_fixSm_ss.gr.in WUSmap_2014_fixSm_ss_M8.in; do
          logf=`echo $inputf | sed 's/in$/txt/' | awk '{print "../logs/log_"$1}'`
          echo "Running $inputf..."
          echo "logfile: $logf"
          $GRID $D_WUS_G/$inputf > $logf
        done

        # WUS-compressional, complement Puget region
	for inputf in noPuget_2014_adSm_r.ch.in noPuget_2014_adSm_r.gr.in noPuget_2014_adSm_r_M8.in noPuget_2014_adSm_ss.ch.in noPuget_2014_adSm_ss.gr.in noPuget_2014_adSm_ss_M8.in noPuget_2014_fixSm_r.ch.in noPuget_2014_fixSm_r.gr.in noPuget_2014_fixSm_r_M8.in noPuget_2014_fixSm_ss.ch.in noPuget_2014_fixSm_ss.gr.in noPuget_2014_fixSm_ss_M8.in puget_2014_r.ch.in puget_2014_r.gr.in puget_2014_r_M8.in puget_2014_ss.ch.in puget_2014_ss.gr.in puget_2014_ss_M8.in; do
          logf=`echo $inputf | sed 's/in$/txt/' | awk '{print "../logs/log_"$1}'`
          echo "Running $inputf..."
          echo "logfile: $logf"
          $GRID $D_WUS_G/$inputf > $logf
        done

        # gridded fault hazard from Bird geodetic model
	for inputf in WUS_zones_PB.in; do
          logf=`echo $inputf | sed 's/in$/txt/' | awk '{print "../logs/log_"$1}'`
          echo "Running $inputf..."
          echo "logfile: $logf"
          $GRID $D_WUS_G/$inputf > $logf
        done

        # WUS shear zones
	for inputf in shear2_2014.in shear3_2014.in shear4_2014.in; do
          logf=`echo $inputf | sed 's/in$/txt/' | awk '{print "../logs/log_"$1}'`
          echo "Running $inputf..."
          echo "logfile: $logf"
          $GRID $D_WUS_G/$inputf > $logf
        done

        # deep seismicity - CA, coastal OR, puget sound
	for inputf in CAdeep.2014.in CAdeepMmax75.2014.in CAdeepMmax8.2014.in coastalOR_deep.in coastalOR_deep_Mmax75.in coastalOR_deep_Mmax8.in pacnwdeep.2014.in pacnwdeep_Mmax75.2014.in pacnwdeep_Mmax8.2014.in; do
          logf=`echo $inputf | sed 's/in$/txt/' | awk '{print "../logs/log_"$1}'`
          echo "Running $inputf..."
          echo "logfile: $logf"
          $GRID $D_WUS_G/$inputf > $logf
        done

fi

echo "$0 script completed"

