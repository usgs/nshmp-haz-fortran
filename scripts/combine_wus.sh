#!/bin/bash

BIN=../bin
CONF=../conf/combine
OUT=../out/combine

if [ ! -d $OUT ]; then
  mkdir $OUT
fi

COMBO=$BIN/hazallXLv2
REGRID=$BIN/hazinterpnga

# Combine Cascadia source results
$COMBO "$CONF/casc.f.pga"
$COMBO "$CONF/casc.f.1hz"
$COMBO "$CONF/casc.f.5hz"

$REGRID < "$CONF/casc.f.resample"
rm $OUT/tmp.casc.*


# Combine gridded source results
$COMBO "$CONF/wus.g.pga"
$COMBO "$CONF/wus.g.pga.m"
$COMBO "$CONF/wus.g.pga.p"
$COMBO "$CONF/wus.g.1hz"
$COMBO "$CONF/wus.g.1hz.m"
$COMBO "$CONF/wus.g.1hz.p"
$COMBO "$CONF/wus.g.5hz"
$COMBO "$CONF/wus.g.5hz.m"
$COMBO "$CONF/wus.g.5hz.p"

$COMBO "$CONF/wus.g.pga.epimerge"
$COMBO "$CONF/wus.g.1hz.epimerge"
$COMBO "$CONF/wus.g.5hz.epimerge"

$REGRID < "$CONF/wus.g.resample"
rm $OUT/tmp.wus.g.*


# Combine faults source results
$COMBO "$CONF/wus.f.pga"
$COMBO "$CONF/wus.f.pga.m"
$COMBO "$CONF/wus.f.pga.p"
$COMBO "$CONF/wus.f.1hz"
$COMBO "$CONF/wus.f.1hz.m"
$COMBO "$CONF/wus.f.1hz.p"
$COMBO "$CONF/wus.f.5hz"
$COMBO "$CONF/wus.f.5hz.m"
$COMBO "$CONF/wus.f.5hz.p"

$COMBO "$CONF/wus.f.pga.epimerge"
$COMBO "$CONF/wus.f.1hz.epimerge"
$COMBO "$CONF/wus.f.5hz.epimerge"

rm $OUT/tmp.wus.f.*


# Combine gridded and fault sources
$COMBO "$CONF/wus.pga"
$COMBO "$CONF/wus.1hz"
$COMBO "$CONF/wus.5hz"
