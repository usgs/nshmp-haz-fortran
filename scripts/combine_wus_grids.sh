#!/bin/bash
# script info
# 
 
# EXTmap
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_ad.g.1hz
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_ad.g.5hz
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_ad.g.pga
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_fix.g.1hz
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_fix.g.5hz
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_fix.g.pga
# WUSmap
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_ad.g.1hz
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_ad.g.5hz
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_ad.g.pga
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_fix.g.1hz
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_fix.g.5hz
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_fix.g.pga
# noPuget
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_ad.g.1hz
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_ad.g.5hz
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_ad.g.pga
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_fix.g.1hz
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_fix.g.5hz
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_fix.g.pga
# combine fixed, adaptive and noPuget
../bin/hazallXL.v2 ../conf/combine/wus_2014.g.5hz
../bin/hazallXL.v2 ../conf/combine/wus_2014.g.1hz
../bin/hazallXL.v2 ../conf/combine/wus_2014.g.pga
# interpolate
../bin/hazinterpnga < ../conf/combine/wus_2014.g.resample
