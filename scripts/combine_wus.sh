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

# CASCADIA 
# cascadia source results
$COMBO "$CONF/casc.f.pga"
$COMBO "$CONF/casc.f.1hz"
$COMBO "$CONF/casc.f.5hz"
# re-grid wus files
$REGRID < "$CONF/casc.f.resample"
rm $OUT/tmp.casc.*

# WUS
# wus gridded source results
$COMBO "$CONF/wus.g.pga"
$COMBO "$CONF/wus.g.pga.m"
$COMBO "$CONF/wus.g.pga.p"
$COMBO "$CONF/wus.g.1hz"
$COMBO "$CONF/wus.g.1hz.m"
$COMBO "$CONF/wus.g.1hz.p"
$COMBO "$CONF/wus.g.5hz"
$COMBO "$CONF/wus.g.5hz.m"
$COMBO "$CONF/wus.g.5hz.p"
# merge p-/m-/central-branches
$COMBO "$CONF/wus.g.pga.epimerge"
$COMBO "$CONF/wus.g.1hz.epimerge"
$COMBO "$CONF/wus.g.5hz.epimerge"
# re-grid wus files
$REGRID < "$CONF/wus.g.resample"
rm $OUT/tmp.wus.g.*

# WUS FAULTS SOURCE RESULTS
$COMBO "$CONF/wus.f.pga"
$COMBO "$CONF/wus.f.pga.m"
$COMBO "$CONF/wus.f.pga.p"
$COMBO "$CONF/wus.f.1hz"
$COMBO "$CONF/wus.f.1hz.m"
$COMBO "$CONF/wus.f.1hz.p"
$COMBO "$CONF/wus.f.5hz"
$COMBO "$CONF/wus.f.5hz.m"
$COMBO "$CONF/wus.f.5hz.p"
# merge p-/m-/central-branches
$COMBO "$CONF/wus.f.pga.epimerge"
$COMBO "$CONF/wus.f.1hz.epimerge"
$COMBO "$CONF/wus.f.5hz.epimerge"
# remove temporary files
rm $OUT/tmp.wus.f.*
# combine wus gridded and fault sources
$COMBO "$CONF/wus.pga"
$COMBO "$CONF/wus.1hz"
$COMBO "$CONF/wus.5hz"
