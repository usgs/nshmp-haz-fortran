#!/bin/bash
# script info
# 
 
# CHECK INPUT PARAMETERS
reqparams=2
if [ $# -ne $reqparams ]; then
  echo "USAGE: $0 [hazard directory] [GMPE name]"
  echo "(e.g., /Users/mmoschetti/Documents/work/NSHM_calculations/NSHMP/trunk)"
  echo "GMPE name from: cy07, cb-08, ba07"
  exit
fi
based=$1
gmpenm=$2
time1=`date`
echo "time1: $time1"

# set run logicals
# make modified files
log_gen_files=1
log_check_files=1
# hazard calcs
log_run_haz=0
log_run_combine=0
log_run_haz_check=0
# plotting
log_plot=0
log_mail=0

# set vars
scripts_dir=$based/scripts
conf_dir=$based/conf
bin_dir=$based/bin
log_dir=$based/logs/${gmpenm}_logs
ck_dir=$based/dir_ck_files
combine_dir=$based/conf/combine
gridded_dir=$based/conf/WUS/gridded
faults_dir=$based/conf/WUS/faults
faults_dirCA=$based/conf/CA/faults
gridded_dir_rel=../conf/WUS/gridded
faults_dir_rel=../conf/WUS/faults
faults_dir_relCA=../conf/CA/faults
combine_dir_rel=../conf/combine

# make directory for check of all modified files
if [ ! -d $ck_dir ]; then
  mkdir $ck_dir
fi
if [ ! -d $log_dir ]; then
  mkdir $log_dir
fi
 
# MAKE MODIFIED HAZARD FILES 
if (( $log_gen_files )); then
  echo "Prep modified hazard files for GMPE - $gmpenm:"
# fault files
  echo "  fault files..."
# WUS
  cd $faults_dir
  desc=`grep $gmpenm list_gmpes_faults_WUS.txt`
  for tempf in *_temp*.in; do
    outf=`echo $tempf | sed 's/temp/'$gmpenm'/'`
    cat $tempf | sed "s/_gmpe/_$gmpenm/" | sed "s/GMPE/$desc/" > $outf
  done
  if (( $log_check_files )); then
    cp *${gmpenm}* $ck_dir/
  fi
# gridded files - no gridded files, at this time
#  echo "  gridded files..."
#  cd $gridded_dir
# use different GMPE flages for mblg->Mw conversions
#  desc=`grep $gmpenm list_gmpes_gridded.txt`
#  Mw sources
#  desc=`grep $gmpenm list_gmpes_gridded_Mw.txt`
#  for tempf in CEUSchar*_temp.in; do
#    outf=`echo $tempf | sed 's/temp/'${gmpenm}'/'`
#    cat $tempf | sed "s/gmpe/$gmpenm/" | sed "s/GMPE/$desc/" > $outf
#  done
#
#  if (( $log_check_files )); then
#    cp *${gmpenm}* $ck_dir/
#  fi
  
# modify combine files to get hazard sensitivity from GMPEs
#  echo "  combine files..."
#  cd $combine_dir
# combine fault/gridded files
#  for extv in f.pga f.1hz f.5hz; do
#    freqv=`echo $extv | cut -d. -f2`
#    cat wus_gmpe.${extv} | sed 's/gmpe/'${gmpenm}'/' > wus_${gmpenm}.$extv
#    cat wus_${gmpenm}.$extv | sed 's/'${freqv}'/'${freqv}'.m/' > wus_${gmpenm}.${extv}.m
#    cat wus_${gmpenm}.$extv | sed 's/'${freqv}'/'${freqv}'.p/' > wus_${gmpenm}.${extv}.p
#    cat wus.${extv}.epimerge | sed 's/wus/wus_'${gmpenm}'/' > wus_${gmpenm}.${extv}.epimerge 
#    cat wus.${extv}.epimerge | sed 's/\.'${freqv}'/_'${gmpenm}'\.'${freqv}'/' > wus_${gmpenm}.${extv}.epimerge 
#    cat wus.${freqv} | sed 's/wus\.f\.'${freqv}'/wus\.f_'${gmpenm}'\.'${freqv}'/' | sed 's/wus\.'${freqv}'/wus_'${gmpenm}'\.'${freqv}'/' > wus_${gmpenm}.${freqv}
#  done
#  cat wus.pga | sed 's/wus.f/wus_'${gmpenm}'.f/' | sed 's/wus.pga/wus_'${gmpenm}'/' > wus_${gmpenm}.pga
#  cat wus.1hz | sed 's/wus.f/wus_'${gmpenm}'.f/' | sed 's/wus.1hz/wus_'${gmpenm}'/' > wus_${gmpenm}.1hz
#  cat wus.5hz | sed 's/wus.f/wus_'${gmpenm}'.f/' | sed 's/wus.5hz/wus_'${gmpenm}'/' > wus_${gmpenm}.5hz
# map files
#  cat map_ceus_wus_2pc50.pga | sed 's/wus/wus_'${gmpenm}'/' | sed 's/hazard/hazard_'${gmpenm}'/' > map_ceus_wus_2pc50_${gmpenm}.pga 
#  cat map_ceus_wus_2pc50.1hz | sed 's/wus/wus_'${gmpenm}'/' | sed 's/hazard/hazard_'${gmpenm}'/' > map_ceus_wus_2pc50_${gmpenm}.1hz 
#  cat map_ceus_wus_2pc50.5hz | sed 's/wus/wus_'${gmpenm}'/' | sed 's/hazard/hazard_'${gmpenm}'/' > map_ceus_wus_2pc50_${gmpenm}.5hz 
 # 
  if (( $log_check_files )); then
    cp *${gmpenm}* $ck_dir/
  fi
  echo "Completed making all modified hazard files."
fi

# RUN MODIFIED HAZARD FILES 
if (( $log_run_haz || $log_run_haz_check )); then
# fault files
  echo "Starting hazard runs for modified files:"
  cd $scripts_dir
#  echo "scripts directory: $scripts_dir"
  echo "  fault files..."
#WUS
  for modf in `ls $faults_dir | grep _${gmpenm}`; do
    if (( $log_run_haz_check )); then
      echo $bin_dir/hazFXnga13l $faults_dir_rel/$modf 
    fi
    if (( $log_run_haz )); then
      logf=`echo $modf | sed 's/\.in//'`
      $bin_dir/hazFXnga13l $faults_dir_rel/$modf >& $log_dir/${logf}.log
    fi
  done
#CA
  for modf in `ls $faults_dirCA | grep _${gmpenm}`; do
    if (( $log_run_haz_check )); then
      echo $bin_dir/hazFXnga13l $faults_dir_relCA/$modf 
    fi
    if (( $log_run_haz )); then
      logf=`echo $modf | sed 's/\.in//'`
      $bin_dir/hazFXnga13l $faults_dir_relCA/$modf >& $log_dir/${logf}.log
    fi
  done
# gridded files
  echo "  no gridded files to run, at this time..."
#  for modf in `ls $gridded_dir | grep _${gmpenm}`; do
#    if (( $log_run_haz_check )); then
#      echo $bin_dir/hazgridXnga5 $gridded_dir_rel/$modf 
#    fi
#    if (( $log_run_haz )); then
#      logf=`echo $modf | sed 's/\.in//'`
#      $bin_dir/hazgridXnga5 $gridded_dir_rel/$modf >& $log_dir/${logf}.log
#    fi
#  done
fi
time1a=`date`
echo "time1a: $time1 - completed all hazard calcs"

# COMBINING FILES
if (( $log_run_combine || $log_run_haz_check )); then
# combine ceus source results
# fault/gridded sources
  echo "Combining fault and gridded results (WUS/CEUS)."
  if (( $log_run_haz_check )); then
    echo $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.f.pga
    echo $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.f.pga.m
    echo $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.f.pga.p
    echo $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.f.1hz
    echo $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.f.1hz.m
    echo $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.f.1hz.p
    echo $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.f.5hz
    echo $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.f.5hz.m
    echo $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.f.5hz.p
  fi
  if (( $log_run_combine )); then
    $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.f.pga
    $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.f.pga.m
    $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.f.pga.p
    $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.f.1hz
    $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.f.1hz.m
    $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.f.1hz.p
    $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.f.5hz
    $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.f.5hz.m
    $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.f.5hz.p
  fi
# epimerge
  echo "Epimerge and combine gridded/fault files."
  if (( $log_run_haz_check )); then
    echo $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.f.pga.epimerge
    echo $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.f.1hz.epimerge
    echo $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.f.5hz.epimerge
    echo $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.pga
    echo $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.1hz
    echo $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.5hz
  fi
  if (( $log_run_combine )); then
    $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.f.pga.epimerge
    $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.f.1hz.epimerge
    $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.f.5hz.epimerge
    $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.pga
    $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.1hz
    $bin_dir/hazallXLv2 $combine_dir_rel/wus_${gmpenm}.5hz
  fi
# hazard maps
  echo "Extract hazard maps at PGA, 1 and 5 Hz response-spectral acceleration"
  if (( $log_run_haz_check )); then
    echo $bin_dir/hazallXL.v4 $combine_dir_rel/map_ceus_wus_2pc50_${gmpenm}.pga
    echo $bin_dir/hazallXL.v4 $combine_dir_rel/map_ceus_wus_2pc50_${gmpenm}.1hz
    echo $bin_dir/hazallXL.v4 $combine_dir_rel/map_ceus_wus_2pc50_${gmpenm}.5hz
  fi
  if (( $log_run_combine )); then
    $bin_dir/hazallXL.v4 $combine_dir_rel/map_ceus_wus_2pc50_${gmpenm}.pga
    $bin_dir/hazallXL.v4 $combine_dir_rel/map_ceus_wus_2pc50_${gmpenm}.1hz
    $bin_dir/hazallXL.v4 $combine_dir_rel/map_ceus_wus_2pc50_${gmpenm}.5hz
  fi
  echo "Completed all processing. Log files written to directory, $log_dir"
fi

# PLOTS
if (( $log_plot )); then
cd $based
cp misc/haz_AbsDiffPosNeg2.cpt .
cp misc/haz_Ratio_1.cpt .
minlon=-125
maxlon=-105
minlat=24.6
maxlat=50.0
plot_sz=7.5
cpt_diff=haz_AbsDiffPosNeg2.cpt
cpt_ratio=haz_Ratio_1.cpt
# PGA
echo "NSHM08 -/ ${gmpenm}-GMPE, PGA, 2%PE-50y" > text.tmp
textf=text.tmp
haz_map_2008=out/combine/us_hazard.pga.2pc50
haz_map_gmpe=out/combine/us_hazard_${gmpenm}.pga.2pc50
cpt_files="$cpt_diff $cpt_ratio"
coords="$minlon $maxlon $minlat $maxlat"
maps_f="$haz_map_2008 $haz_map_gmpe"
#echo plot_haz_diffRatio_boundaries.gmt $maps_f $textf $coords $cpt_files $plot_sz
#plot_haz_diffRatio_boundaries.gmt $maps_f $textf $coords $cpt_files $plot_sz
echo plot_haz_diffRatio_boundaries2.gmt $maps_f $textf $coords $cpt_files $plot_sz
plot_haz_diffRatio_boundaries2.gmt $maps_f $textf $coords $cpt_files $plot_sz
mv pl_hazard_diff_ratio.ps gmpe_tests_wus/pl_haz_${gmpenm}_pga_2pc50.ps
# 1 Hz
echo "NSHM08 -/ ${gmpenm}-GMPE, 1-Hz, 2%PE-50y" > text.tmp
textf=text.tmp
haz_map_2008=out/combine/us_hazard.1hz.2pc50
haz_map_gmpe=out/combine/us_hazard_${gmpenm}.1hz.2pc50
cpt_files="$cpt_diff $cpt_ratio"
coords="$minlon $maxlon $minlat $maxlat"
maps_f="$haz_map_2008 $haz_map_gmpe"
echo plot_haz_diffRatio_boundaries.gmt $maps_f $textf $coords $cpt_files $plot_sz
#plot_haz_diffRatio_boundaries.gmt $maps_f $textf $coords $cpt_files $plot_sz
plot_haz_diffRatio_boundaries2.gmt $maps_f $textf $coords $cpt_files $plot_sz
mv pl_hazard_diff_ratio.ps gmpe_tests_wus/pl_haz_${gmpenm}_1hz_2pc50.ps
# 5 Hz
echo "NSHM08 -/ ${gmpenm}-GMPE, 5-Hz, 2%PE-50y" > text.tmp
textf=text.tmp
haz_map_2008=out/combine/us_hazard.5hz.2pc50
haz_map_gmpe=out/combine/us_hazard_${gmpenm}.5hz.2pc50
cpt_files="$cpt_diff $cpt_ratio"
coords="$minlon $maxlon $minlat $maxlat"
maps_f="$haz_map_2008 $haz_map_gmpe"
echo plot_haz_diffRatio_boundaries.gmt $maps_f $textf $coords $cpt_files $plot_sz
#plot_haz_diffRatio_boundaries.gmt $maps_f $textf $coords $cpt_files $plot_sz
plot_haz_diffRatio_boundaries2.gmt $maps_f $textf $coords $cpt_files $plot_sz
mv pl_hazard_diff_ratio.ps gmpe_tests_wus/pl_haz_${gmpenm}_5hz_2pc50.ps
fi


if (( $log_mail )); then
time2=`date`

mail -s "`basename $0` - $gmpenm" mmoschetti@usgs.gov << ENDMAIL
$0 script completed

Start time: $time1
End time: $time2

ENDMAIL
fi
