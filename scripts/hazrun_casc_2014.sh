#!/bin/bash

time1=`date`

LOG=../logs
BIN=../bin
CONF=../conf

SUB=$BIN/hazSUBX

CASC_CONF=$CONF/CASC
CASC_SUB=1

# All Subduction
if [ $CASC_SUB == 1 ]; then

        # sub0
        $SUB $CASC_CONF/sub0_GRb0_bot.in > $LOG/log_sub0_GRb0_bot.txt
        $SUB $CASC_CONF/sub0_GRb0_mid.in > $LOG/log_sub0_GRb0_mid.txt
        $SUB $CASC_CONF/sub0_GRb0_top.in > $LOG/log_sub0_GRb0_top.txt
        $SUB $CASC_CONF/sub0_GRb1_bot.in > $LOG/log_sub0_GRb1_bot.txt
        $SUB $CASC_CONF/sub0_GRb1_mid.in > $LOG/log_sub0_GRb1_mid.txt
        $SUB $CASC_CONF/sub0_GRb1_top.in > $LOG/log_sub0_GRb1_top.txt
        $SUB $CASC_CONF/sub0_ch_bot.in > $LOG/log_sub0_ch_bot.txt
        $SUB $CASC_CONF/sub0_ch_mid.in > $LOG/log_sub0_ch_mid.txt
        $SUB $CASC_CONF/sub0_ch_top.in > $LOG/log_sub0_ch_top.txt
        # sub1
        $SUB $CASC_CONF/sub1_GRb0_bot.in > $LOG/log_sub1_GRb0_bot.txt
        $SUB $CASC_CONF/sub1_GRb0_mid.in > $LOG/log_sub1_GRb0_mid.txt
        $SUB $CASC_CONF/sub1_GRb0_top.in > $LOG/log_sub1_GRb0_top.txt
        $SUB $CASC_CONF/sub1_GRb1_bot.in > $LOG/log_sub1_GRb1_bot.txt
        $SUB $CASC_CONF/sub1_GRb1_mid.in > $LOG/log_sub1_GRb1_mid.txt
        $SUB $CASC_CONF/sub1_GRb1_top.in > $LOG/log_sub1_GRb1_top.txt
        $SUB $CASC_CONF/sub1_ch_bot.in > $LOG/log_sub1_ch_bot.txt
        $SUB $CASC_CONF/sub1_ch_mid.in > $LOG/log_sub1_ch_mid.txt
        $SUB $CASC_CONF/sub1_ch_top.in > $LOG/log_sub1_ch_top.txt
        # sub2
        $SUB $CASC_CONF/sub2_ch_bot.in > $LOG/log_sub2_ch_bot.txt
        $SUB $CASC_CONF/sub2_ch_mid.in > $LOG/log_sub2_ch_mid.txt
        $SUB $CASC_CONF/sub2_ch_top.in > $LOG/log_sub2_ch_top.txt
        # sub3
        $SUB $CASC_CONF/sub3_ch_bot.in > $LOG/log_sub3_ch_bot.txt
        $SUB $CASC_CONF/sub3_ch_mid.in > $LOG/log_sub3_ch_mid.txt
        $SUB $CASC_CONF/sub3_ch_top.in > $LOG/log_sub3_ch_top.txt
        # sub4
        $SUB $CASC_CONF/sub4_ch_bot.in > $LOG/log_sub4_ch_bot.txt
        $SUB $CASC_CONF/sub4_ch_mid.in > $LOG/log_sub4_ch_mid.txt
        $SUB $CASC_CONF/sub4_ch_top.in > $LOG/log_sub4_ch_top.txt
fi

