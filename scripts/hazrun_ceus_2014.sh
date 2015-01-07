#!/bin/bash

time1=`date`
LOG=../logs
BIN=../bin
CONF=../conf

FAULT=$BIN/hazFXnga13l
GRID=$BIN/hazgridXnga13l

CEUS_CONF=$CONF/CEUS
CEUS_FAULTS=1
CEUS_BACKGRND=1

# CEUS Faults
if [ $CEUS_FAULTS == 1 ]; then
	D_CEUS_F=$CEUS_CONF/faults

	# Cheraw, Meers faults
        $FAULT $D_CEUS_F/CEUScm-meers_2014.in > $LOG/log_CEUScm-meers_2014.txt
        $FAULT $D_CEUS_F/CEUScm-recur_2014.in > $LOG/log_CEUScm-recur_2014.txt
        $FAULT $D_CEUS_F/CEUScm-srchar_2014.in > $LOG/log_CEUScm-srchar_2014.txt
        $FAULT $D_CEUS_F/CEUScm-srgr_2014.in > $LOG/log_CEUScm-srgr_2014.txt
        $FAULT $D_CEUS_F/CEUScm2014.in > $LOG/log_CEUScm2014.txt

	# NMSZ, CEUS-SSC-N
        $FAULT $D_CEUS_F/NMFS_RFT.RLME.in > $LOG/log_NMFS_RFT.RLME.txt
        $FAULT $D_CEUS_F/NMFS_RLME_clu.in > $LOG/log_NMFS_RLME_clu.txt
	# NMSZ, un-clustered
        $FAULT $D_CEUS_F/NMSZnocl.1000.2014.in > $LOG/log_NMSZnocl.1000.2014.txt
        $FAULT $D_CEUS_F/NMSZnocl.500.2014.in > $LOG/log_NMSZnocl.500.2014.txt
        $FAULT $D_CEUS_F/NMSZnocl.50000.2014.in > $LOG/log_NMSZnocl.50000.2014.txt
	# NMSZ, clustered
        $FAULT $D_CEUS_F/newmad2014.1000.cluster.in > $LOG/log_newmad2014.1000.cluster.txt
        $FAULT $D_CEUS_F/newmad2014.1500.cluster.in > $LOG/log_newmad2014.1500.cluster.txt
        $FAULT $D_CEUS_F/newmad2014.500.cluster.0p3.in > $LOG/log_newmad2014.500.cluster.0p3.txt
        $FAULT $D_CEUS_F/newmad2014.500.cluster.in > $LOG/log_newmad2014.500.cluster.txt
        $FAULT $D_CEUS_F/newmad2014.50000.cluster.in > $LOG/log_newmad2014.50000.cluster.txt
        $FAULT $D_CEUS_F/newmad2014.750.cluster.in > $LOG/log_newmad2014.750.cluster.txt
fi

# CEUS Background 
if [ $CEUS_BACKGRND == 1 ]; then
	D_CEUS_G=$CEUS_CONF/gridded

	# gridded-seismicity
        $GRID $D_CEUS_G/CEUS_adaptGridded_2014_2zone.in > $LOG/log_CEUS_adaptGridded_2014_2zone.txt
        $GRID $D_CEUS_G/CEUS_adaptGridded_2014_2zone_0p3.in > $LOG/log_CEUS_adaptGridded_2014_2zone_0p3.txt
        $GRID $D_CEUS_G/CEUS_adaptGridded_2014_4zone.in > $LOG/log_CEUS_adaptGridded_2014_4zone.txt
        $GRID $D_CEUS_G/CEUS_fixRGridded_2014_2zone.in > $LOG/log_CEUS_fixRGridded_2014_2zone.txt
        $GRID $D_CEUS_G/CEUS_fixRGridded_2014_4zone.in > $LOG/log_CEUS_fixRGridded_2014_4zone.txt

	# RLMEs
	# Charleston
        $GRID $D_CEUS_G/CEUSchar_2014_l.ssc67.in > $LOG/log_CEUSchar_2014_l.ssc67.txt
        $GRID $D_CEUS_G/CEUSchar_2014_l.ssc69.in > $LOG/log_CEUSchar_2014_l.ssc69.txt
        $GRID $D_CEUS_G/CEUSchar_2014_l.ssc71.in > $LOG/log_CEUSchar_2014_l.ssc71.txt
        $GRID $D_CEUS_G/CEUSchar_2014_l.ssc73.in > $LOG/log_CEUSchar_2014_l.ssc73.txt
        $GRID $D_CEUS_G/CEUSchar_2014_l.ssc75.in > $LOG/log_CEUSchar_2014_l.ssc75.txt
        $GRID $D_CEUS_G/CEUSchar_2014_n.ssc67.in > $LOG/log_CEUSchar_2014_n.ssc67.txt
        $GRID $D_CEUS_G/CEUSchar_2014_n.ssc69.in > $LOG/log_CEUSchar_2014_n.ssc69.txt
        $GRID $D_CEUS_G/CEUSchar_2014_n.ssc71.in > $LOG/log_CEUSchar_2014_n.ssc71.txt
        $GRID $D_CEUS_G/CEUSchar_2014_n.ssc73.in > $LOG/log_CEUSchar_2014_n.ssc73.txt
        $GRID $D_CEUS_G/CEUSchar_2014_n.ssc75.in > $LOG/log_CEUSchar_2014_n.ssc75.txt
        $GRID $D_CEUS_G/CEUSchar_2014_r.ssc67.in > $LOG/log_CEUSchar_2014_r.ssc67.txt
        $GRID $D_CEUS_G/CEUSchar_2014_r.ssc69.in > $LOG/log_CEUSchar_2014_r.ssc69.txt
        $GRID $D_CEUS_G/CEUSchar_2014_r.ssc71.in > $LOG/log_CEUSchar_2014_r.ssc71.txt
        $GRID $D_CEUS_G/CEUSchar_2014_r.ssc73.in > $LOG/log_CEUSchar_2014_r.ssc73.txt
        $GRID $D_CEUS_G/CEUSchar_2014_r.ssc75.in > $LOG/log_CEUSchar_2014_r.ssc75.txt
	# Charlevoix
        $GRID $D_CEUS_G/Chlvx_RLME_2014.675.in > $LOG/log_Chlvx_RLME_2014.675.txt
        $GRID $D_CEUS_G/Chlvx_RLME_2014.700.in > $LOG/log_Chlvx_RLME_2014.700.txt
        $GRID $D_CEUS_G/Chlvx_RLME_2014.725.in > $LOG/log_Chlvx_RLME_2014.725.txt
        $GRID $D_CEUS_G/Chlvx_RLME_2014.750.in > $LOG/log_Chlvx_RLME_2014.750.txt
	# Commerce Geophysical Lineament
        $GRID $D_CEUS_G/Commerce_RLME_2014.670.in > $LOG/log_Commerce_RLME_2014.670.txt
        $GRID $D_CEUS_G/Commerce_RLME_2014.690.in > $LOG/log_Commerce_RLME_2014.690.txt
        $GRID $D_CEUS_G/Commerce_RLME_2014.710.in > $LOG/log_Commerce_RLME_2014.710.txt
        $GRID $D_CEUS_G/Commerce_RLME_2014.730.in > $LOG/log_Commerce_RLME_2014.730.txt
        $GRID $D_CEUS_G/Commerce_RLME_2014.750.in > $LOG/log_Commerce_RLME_2014.750.txt
        $GRID $D_CEUS_G/Commerce_RLME_2014.770.in > $LOG/log_Commerce_RLME_2014.770.txt
	# East Rift Margin
        $GRID $D_CEUS_G/ERM-N_RLME_2014.670.in > $LOG/log_ERM-N_RLME_2014.670.txt
        $GRID $D_CEUS_G/ERM-N_RLME_2014.690.in > $LOG/log_ERM-N_RLME_2014.690.txt
        $GRID $D_CEUS_G/ERM-N_RLME_2014.710.in > $LOG/log_ERM-N_RLME_2014.710.txt
        $GRID $D_CEUS_G/ERM-N_RLME_2014.740.in > $LOG/log_ERM-N_RLME_2014.740.txt
        $GRID $D_CEUS_G/ERM-S1_RLME_2014.670.in > $LOG/log_ERM-S1_RLME_2014.670.txt
        $GRID $D_CEUS_G/ERM-S1_RLME_2014.690.in > $LOG/log_ERM-S1_RLME_2014.690.txt
        $GRID $D_CEUS_G/ERM-S1_RLME_2014.710.in > $LOG/log_ERM-S1_RLME_2014.710.txt
        $GRID $D_CEUS_G/ERM-S1_RLME_2014.730.in > $LOG/log_ERM-S1_RLME_2014.730.txt
        $GRID $D_CEUS_G/ERM-S1_RLME_2014.750.in > $LOG/log_ERM-S1_RLME_2014.750.txt
        $GRID $D_CEUS_G/ERM-S1_RLME_2014.770.in > $LOG/log_ERM-S1_RLME_2014.770.txt
        $GRID $D_CEUS_G/ERM-S2_RLME_2014.670.in > $LOG/log_ERM-S2_RLME_2014.670.txt
        $GRID $D_CEUS_G/ERM-S2_RLME_2014.690.in > $LOG/log_ERM-S2_RLME_2014.690.txt
        $GRID $D_CEUS_G/ERM-S2_RLME_2014.710.in > $LOG/log_ERM-S2_RLME_2014.710.txt
        $GRID $D_CEUS_G/ERM-S2_RLME_2014.730.in > $LOG/log_ERM-S2_RLME_2014.730.txt
        $GRID $D_CEUS_G/ERM-S2_RLME_2014.750.in > $LOG/log_ERM-S2_RLME_2014.750.txt
        $GRID $D_CEUS_G/ERM-S2_RLME_2014.770.in > $LOG/log_ERM-S2_RLME_2014.770.txt
	# Marianna
        $GRID $D_CEUS_G/Marianna_RLME_2014.670.in > $LOG/log_Marianna_RLME_2014.670.txt
        $GRID $D_CEUS_G/Marianna_RLME_2014.690.in > $LOG/log_Marianna_RLME_2014.690.txt
        $GRID $D_CEUS_G/Marianna_RLME_2014.710.in > $LOG/log_Marianna_RLME_2014.710.txt
        $GRID $D_CEUS_G/Marianna_RLME_2014.730.in > $LOG/log_Marianna_RLME_2014.730.txt
        $GRID $D_CEUS_G/Marianna_RLME_2014.750.in > $LOG/log_Marianna_RLME_2014.750.txt
        $GRID $D_CEUS_G/Marianna_RLME_2014.770.in > $LOG/log_Marianna_RLME_2014.770.txt
	# Wabash
        $GRID $D_CEUS_G/Wabash_RLME_2014.675.in > $LOG/log_Wabash_RLME_2014.675.txt
        $GRID $D_CEUS_G/Wabash_RLME_2014.700.in > $LOG/log_Wabash_RLME_2014.700.txt
        $GRID $D_CEUS_G/Wabash_RLME_2014.725.in > $LOG/log_Wabash_RLME_2014.725.txt
        $GRID $D_CEUS_G/Wabash_RLME_2014.750.in > $LOG/log_Wabash_RLME_2014.750.txt
fi
