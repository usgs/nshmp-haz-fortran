#!/bin/bash

BIN=../bin
CONF=../conf/combine
OUT=../out/combine

if [ ! -d $OUT ]; then
  mkdir $OUT
fi

COMBO=$BIN/hazallXLv2
COMBO_MAP=$BIN/hazallXL.v4
REGRID=$BIN/hazinterpnga

# CEUS
# ceus source results
$COMBO "$CONF/ceus.f.pga"
$COMBO "$CONF/ceus.f.1hz"
$COMBO "$CONF/ceus.f.5hz"
# re-sample CEUS
$REGRID < "$CONF/ceus.all.resample"
rm $OUT/tmp.ceus.*

# combine the WUS and CEUS combined files
# PGA 2% in 50 yrs
$COMBO_MAP "$CONF/map_ceus_wus_2pc50.pga"
$COMBO_MAP "$CONF/map_ceus_wus_2pc50.1hz"
$COMBO_MAP "$CONF/map_ceus_wus_2pc50.5hz"
