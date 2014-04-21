#!/bin/bash
# script info
# 

run_WUS_faults=0
run_WUS_grids=0
run_CEUS_faults=1
run_CA=0
run_CASC=0
run_all_hazard_curves=1
 
# WUS faults
if (( $run_WUS_faults )); then
# central - branch
../bin/hazallXL.v2 ../conf/combine/combine_Geo.5hz
../bin/hazallXL.v2 ../conf/combine/combine_Geo.1hz
../bin/hazallXL.v2 ../conf/combine/combine_Geo.pga
../bin/hazallXL.v2 ../conf/combine/combine_Zeng.5hz
../bin/hazallXL.v2 ../conf/combine/combine_Zeng.1hz
../bin/hazallXL.v2 ../conf/combine/combine_Zeng.pga
../bin/hazallXL.v2 ../conf/combine/combine_Bird.5hz
../bin/hazallXL.v2 ../conf/combine/combine_Bird.1hz
../bin/hazallXL.v2 ../conf/combine/combine_Bird.pga
# combine geologic and geodetic
../bin/hazallXL.v2 ../conf/combine/wus_2014.f.5hz 
../bin/hazallXL.v2 ../conf/combine/wus_2014.f.1hz 
../bin/hazallXL.v2 ../conf/combine/wus_2014.f.pga
# plus - branch
../bin/hazallXL.v2 ../conf/combine/combine_Geo.5hz.p
../bin/hazallXL.v2 ../conf/combine/combine_Geo.1hz.p
../bin/hazallXL.v2 ../conf/combine/combine_Geo.pga.p
../bin/hazallXL.v2 ../conf/combine/combine_Zeng.5hz.p
../bin/hazallXL.v2 ../conf/combine/combine_Zeng.1hz.p
../bin/hazallXL.v2 ../conf/combine/combine_Zeng.pga.p
../bin/hazallXL.v2 ../conf/combine/combine_Bird.5hz.p
../bin/hazallXL.v2 ../conf/combine/combine_Bird.1hz.p
../bin/hazallXL.v2 ../conf/combine/combine_Bird.pga.p
# combine geologic and geodetic
../bin/hazallXL.v2 ../conf/combine/wus_2014.f.5hz.p
../bin/hazallXL.v2 ../conf/combine/wus_2014.f.1hz.p
../bin/hazallXL.v2 ../conf/combine/wus_2014.f.pga.p
# negative - branch
../bin/hazallXL.v2 ../conf/combine/combine_Geo.5hz.m
../bin/hazallXL.v2 ../conf/combine/combine_Geo.1hz.m
../bin/hazallXL.v2 ../conf/combine/combine_Geo.pga.m
../bin/hazallXL.v2 ../conf/combine/combine_Zeng.5hz.m
../bin/hazallXL.v2 ../conf/combine/combine_Zeng.1hz.m
../bin/hazallXL.v2 ../conf/combine/combine_Zeng.pga.m
../bin/hazallXL.v2 ../conf/combine/combine_Bird.5hz.m
../bin/hazallXL.v2 ../conf/combine/combine_Bird.1hz.m
../bin/hazallXL.v2 ../conf/combine/combine_Bird.pga.m
# combine geologic and geodetic
../bin/hazallXL.v2 ../conf/combine/wus_2014.f.5hz.m
../bin/hazallXL.v2 ../conf/combine/wus_2014.f.1hz.m
../bin/hazallXL.v2 ../conf/combine/wus_2014.f.pga.m

# merge branches
../bin/hazallXL.v2 ../conf/combine/wus_2014.f.5hz.epimerge
../bin/hazallXL.v2 ../conf/combine/wus_2014.f.1hz.epimerge
../bin/hazallXL.v2 ../conf/combine/wus_2014.f.pga.epimerge
fi

# WUS grids
if (( $run_WUS_grids )); then
# options for addnl epistemic uncertainty branches
central_branch=1
m_branch=1
p_branch=1
# EXTmap
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_ad.g.1hz
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_ad.g.5hz
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_ad.g.pga
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_fix.g.1hz
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_fix.g.5hz
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_fix.g.pga
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_ad_M8.g.1hz
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_ad_M8.g.5hz
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_ad_M8.g.pga
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_fix_M8.g.1hz
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_fix_M8.g.5hz
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_fix_M8.g.pga
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_ad.g.1hz.p
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_ad.g.5hz.p
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_ad.g.pga.p
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_fix.g.1hz.p
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_fix.g.5hz.p
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_fix.g.pga.p
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_ad_M8.g.1hz.p
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_ad_M8.g.5hz.p
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_ad_M8.g.pga.p
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_fix_M8.g.1hz.p
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_fix_M8.g.5hz.p
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_fix_M8.g.pga.p
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_ad.g.1hz.m
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_ad.g.5hz.m
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_ad.g.pga.m
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_fix.g.1hz.m
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_fix.g.5hz.m
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_fix.g.pga.m
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_ad_M8.g.1hz.m
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_ad_M8.g.5hz.m
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_ad_M8.g.pga.m
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_fix_M8.g.1hz.m
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_fix_M8.g.5hz.m
../bin/hazallXL.v2 ../conf/combine/wus_EXTmap_fix_M8.g.pga.m
# WUSmap
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_ad.g.1hz
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_ad.g.5hz
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_ad.g.pga
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_fix.g.1hz
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_fix.g.5hz
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_fix.g.pga
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_ad_M8.g.1hz
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_ad_M8.g.5hz
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_ad_M8.g.pga
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_fix_M8.g.1hz
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_fix_M8.g.5hz
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_fix_M8.g.pga
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_ad.g.1hz.p
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_ad.g.5hz.p
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_ad.g.pga.p
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_fix.g.1hz.p
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_fix.g.5hz.p
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_fix.g.pga.p
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_ad_M8.g.1hz.p
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_ad_M8.g.5hz.p
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_ad_M8.g.pga.p
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_fix_M8.g.1hz.p
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_fix_M8.g.5hz.p
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_fix_M8.g.pga.p
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_ad.g.1hz.m
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_ad.g.5hz.m
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_ad.g.pga.m
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_fix.g.1hz.m
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_fix.g.5hz.m
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_fix.g.pga.m
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_ad_M8.g.1hz.m
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_ad_M8.g.5hz.m
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_ad_M8.g.pga.m
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_fix_M8.g.1hz.m
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_fix_M8.g.5hz.m
../bin/hazallXL.v2 ../conf/combine/wus_WUSmap_fix_M8.g.pga.m
# noPuget
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_ad.g.1hz
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_ad.g.5hz
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_ad.g.pga
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_fix.g.1hz
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_fix.g.5hz
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_fix.g.pga
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_ad_M8.g.1hz
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_ad_M8.g.5hz
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_ad_M8.g.pga
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_fix_M8.g.1hz
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_fix_M8.g.5hz
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_fix_M8.g.pga
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_ad.g.1hz.p
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_ad.g.5hz.p
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_ad.g.pga.p
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_fix.g.1hz.p
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_fix.g.5hz.p
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_fix.g.pga.p
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_ad_M8.g.1hz.p
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_ad_M8.g.5hz.p
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_ad_M8.g.pga.p
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_fix_M8.g.1hz.p
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_fix_M8.g.5hz.p
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_fix_M8.g.pga.p
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_ad.g.1hz.m
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_ad.g.5hz.m
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_ad.g.pga.m
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_fix.g.1hz.m
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_fix.g.5hz.m
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_fix.g.pga.m
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_ad_M8.g.1hz.m
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_ad_M8.g.5hz.m
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_ad_M8.g.pga.m
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_fix_M8.g.1hz.m
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_fix_M8.g.5hz.m
../bin/hazallXL.v2 ../conf/combine/wus_noPuget_fix_M8.g.pga.m
# puget
../bin/hazallXL.v2 ../conf/combine/wus_puget.g.5hz
../bin/hazallXL.v2 ../conf/combine/wus_puget.g.1hz
../bin/hazallXL.v2 ../conf/combine/wus_puget.g.pga
../bin/hazallXL.v2 ../conf/combine/wus_puget_M8.g.1hz
../bin/hazallXL.v2 ../conf/combine/wus_puget_M8.g.5hz
../bin/hazallXL.v2 ../conf/combine/wus_puget_M8.g.pga
../bin/hazallXL.v2 ../conf/combine/wus_puget.g.5hz.p
../bin/hazallXL.v2 ../conf/combine/wus_puget.g.1hz.p
../bin/hazallXL.v2 ../conf/combine/wus_puget.g.pga.p
../bin/hazallXL.v2 ../conf/combine/wus_puget_M8.g.1hz.p
../bin/hazallXL.v2 ../conf/combine/wus_puget_M8.g.5hz.p
../bin/hazallXL.v2 ../conf/combine/wus_puget_M8.g.pga.p
../bin/hazallXL.v2 ../conf/combine/wus_puget.g.5hz.m
../bin/hazallXL.v2 ../conf/combine/wus_puget.g.1hz.m
../bin/hazallXL.v2 ../conf/combine/wus_puget.g.pga.m
../bin/hazallXL.v2 ../conf/combine/wus_puget_M8.g.1hz.m
../bin/hazallXL.v2 ../conf/combine/wus_puget_M8.g.5hz.m
../bin/hazallXL.v2 ../conf/combine/wus_puget_M8.g.pga.m
# shear zones
../bin/hazallXL.v2 ../conf/combine/shear_zones.g.1hz
../bin/hazallXL.v2 ../conf/combine/shear_zones.g.5hz
../bin/hazallXL.v2 ../conf/combine/shear_zones.g.pga
../bin/hazallXL.v2 ../conf/combine/shear_zones.g.1hz.m
../bin/hazallXL.v2 ../conf/combine/shear_zones.g.5hz.m
../bin/hazallXL.v2 ../conf/combine/shear_zones.g.pga.m
../bin/hazallXL.v2 ../conf/combine/shear_zones.g.1hz.p
../bin/hazallXL.v2 ../conf/combine/shear_zones.g.5hz.p
../bin/hazallXL.v2 ../conf/combine/shear_zones.g.pga.p
# combine all WUS grid sources, except shear zones 
../bin/hazallXL.v2 ../conf/combine/wus_2014.g.5hz
../bin/hazallXL.v2 ../conf/combine/wus_2014.g.1hz
../bin/hazallXL.v2 ../conf/combine/wus_2014.g.pga
../bin/hazallXL.v2 ../conf/combine/wus_2014.g.5hz.m
../bin/hazallXL.v2 ../conf/combine/wus_2014.g.1hz.m
../bin/hazallXL.v2 ../conf/combine/wus_2014.g.pga.m
../bin/hazallXL.v2 ../conf/combine/wus_2014.g.5hz.p
../bin/hazallXL.v2 ../conf/combine/wus_2014.g.1hz.p
../bin/hazallXL.v2 ../conf/combine/wus_2014.g.pga.p

# interpolate 0.1 degree to 0.05 degree
../bin/hazinterpnga < ../conf/combine/wus_shear.all.resample
../bin/hazinterpnga < ../conf/combine/wus_2014.g.resample

# merge branches
../bin/hazallXL.v2 ../conf/combine/wus_2014.g.5hz.epimerge
../bin/hazallXL.v2 ../conf/combine/wus_2014.g.1hz.epimerge
../bin/hazallXL.v2 ../conf/combine/wus_2014.g.pga.epimerge
../bin/hazallXL.v2 ../conf/combine/shear.5hz.epimerge
../bin/hazallXL.v2 ../conf/combine/shear.1hz.epimerge
../bin/hazallXL.v2 ../conf/combine/shear.pga.epimerge
fi

# CA model - from PP
if (( $run_CA )); then
#../bin/hazinterpnga < ../conf/combine/ca.all.resample
echo "Nothing to do for run_CA option"
fi

# CEUS 
if (( $run_CEUS_faults )); then
../bin/hazallXL.v2 ../conf/combine/ceus_Chlvx_2014.f.5hz
../bin/hazallXL.v2 ../conf/combine/ceus_Chlvx_2014.f.1hz
../bin/hazallXL.v2 ../conf/combine/ceus_Chlvx_2014.f.pga
../bin/hazallXL.v2 ../conf/combine/ceus_Charleston_2014.f.5hz
../bin/hazallXL.v2 ../conf/combine/ceus_Charleston_2014.f.1hz
../bin/hazallXL.v2 ../conf/combine/ceus_Charleston_2014.f.pga
../bin/hazallXL.v2 ../conf/combine/ceus_Commerce_2014.f.5hz
../bin/hazallXL.v2 ../conf/combine/ceus_Commerce_2014.f.1hz
../bin/hazallXL.v2 ../conf/combine/ceus_Commerce_2014.f.pga
../bin/hazallXL.v2 ../conf/combine/ceus_ERM-N_2014.f.5hz
../bin/hazallXL.v2 ../conf/combine/ceus_ERM-N_2014.f.1hz
../bin/hazallXL.v2 ../conf/combine/ceus_ERM-N_2014.f.pga
../bin/hazallXL.v2 ../conf/combine/ceus_ERM-S_2014.f.5hz
../bin/hazallXL.v2 ../conf/combine/ceus_ERM-S_2014.f.1hz
../bin/hazallXL.v2 ../conf/combine/ceus_ERM-S_2014.f.pga
../bin/hazallXL.v2 ../conf/combine/ceus_Marianna_2014.f.5hz
../bin/hazallXL.v2 ../conf/combine/ceus_Marianna_2014.f.1hz
../bin/hazallXL.v2 ../conf/combine/ceus_Marianna_2014.f.pga
../bin/hazallXL.v2 ../conf/combine/ceus_Wabash_2014.f.5hz
../bin/hazallXL.v2 ../conf/combine/ceus_Wabash_2014.f.1hz
../bin/hazallXL.v2 ../conf/combine/ceus_Wabash_2014.f.pga
../bin/hazallXL.v2 ../conf/combine/ceus_newmad_clu_2014.f.5hz
../bin/hazallXL.v2 ../conf/combine/ceus_newmad_clu_2014.f.1hz
../bin/hazallXL.v2 ../conf/combine/ceus_newmad_clu_2014.f.pga
../bin/hazallXL.v2 ../conf/combine/ceus_newmad_unclu_2014.f.5hz
../bin/hazallXL.v2 ../conf/combine/ceus_newmad_unclu_2014.f.1hz
../bin/hazallXL.v2 ../conf/combine/ceus_newmad_unclu_2014.f.pga
../bin/hazallXL.v2 ../conf/combine/ceus_NMFS_clu.f.5hz
../bin/hazallXL.v2 ../conf/combine/ceus_NMFS_clu.f.1hz
../bin/hazallXL.v2 ../conf/combine/ceus_NMFS_clu.f.pga
# combine all sources - faults and grids
../bin/hazallXL.v2 ../conf/combine/ceus_2014.f.5hz
../bin/hazallXL.v2 ../conf/combine/ceus_2014.f.1hz
../bin/hazallXL.v2 ../conf/combine/ceus_2014.f.pga
# resample
../bin/hazinterpnga < ../conf/combine/ceus.all.resample
fi

# CASC
if (( $run_CASC )); then
../bin/hazallXL.v2 ../conf/combine/casc_sub2014.1hz
../bin/hazallXL.v2 ../conf/combine/casc_sub2014.5hz
../bin/hazallXL.v2 ../conf/combine/casc_sub2014.pga
../bin/hazinterpnga < ../conf/combine/casc_2014.f.resample
fi

# all-US hazard curves
if (( $run_all_hazard_curves )); then
../bin/hazallXL.v4 ../conf/combine/us_hazard_curves.5hz
../bin/hazallXL.v4 ../conf/combine/us_hazard_curves.1hz
../bin/hazallXL.v4 ../conf/combine/us_hazard_curves.pga
fi


