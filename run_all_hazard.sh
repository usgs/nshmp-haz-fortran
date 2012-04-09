#!/bin/bash
# script to run hazard calculations for lower-48. 
# This repeats the calculations that were used for the 2008, Version III 
# National Seismic Hazard maps. 
# Output files are written to the directories ./out and ./out/combine
# Figures for PGA, 1 and 5 Hz (2% in 50 years) are written to ./figs
# 
 
# CHECK INPUT PARAMETERS
reqparams=2
if [ $# -ne $reqparams ]; then
  echo "USAGE: $0 [path to NSHMP/trunk directory] [email address for final notification]"
  echo "Email will be sent to user when processing is completed"
  exit
fi
trunk_path=$1
email_address=$2

# 
scripts_path=$trunk_path/scripts
time1=`date`

# change to scripts directory
if [ -d $scripts_path ]; then
  cd $scripts_path
else
  echo "$scripts path does not exist"
  exit
fi
 
# create the out and logs directories for output 
if [ ! -d ../out ]; then
  mkdir ../out
fi
if [ ! -d ../logs ]; then
  mkdir ../logs
fi

# re-compile binaries
if [ ! -d ../bin ]; then
  mkdir ../bin
fi
cd ..
make
cd $scripts_path

#check that meanrjb.bin file exists in the bin directory
# option 1 is for WC SRL - used for all gridded seismicity in 2008 maps
cd ../bin
if [ ! -f rbjmean.bin ]; then
getmeanrjf << END
1
END
  mv rjbmean.bin.srl meanrjb.bin
  mv rjbmean.dat.srl meanrjb_srl.dat
fi

# cd to scripts directory for consistency with hazard input files
cd $scripts_path

# run hazard calculation scripts
hazrun_wus.sh >& ../logs/log_hazrun_wus.txt
hazrun_ca.sh >& ../logs/log_hazrun_ca.txt
hazrun_casc.sh >& ../logs/log_hazrun_casc.txt
hazrun_ceus.sh >& ../logs/log_hazrun_ceus.txt

# combine the WUS output files, CEUS output files, merge WUS-CEUS
combine_wus.sh >& ../logs/log_combine_wus.txt
combine_ceus_wus.sh >& ../logs/log_combine_ceus_wus.txt

# plot 2% in 50 yrs results for and convert to pdf
plot_haz_maps_2pc50.gmt ../out/combine/us_hazard.pga.2pc50
plot_haz_maps_2pc50.gmt ../out/combine/us_hazard.1hz.2pc50
plot_haz_maps_2pc50.gmt ../out/combine/us_hazard.5hz.2pc50
if [ ! -d ../figs ]; then
  mkdir ../figs
fi
ps2pdf pl_us_hazard_pga_2pc50.ps ../figs/pl_us_hazard_pga_2pc50.pdf
mv pl_us_hazard_pga_2pc50.ps ../figs/
ps2pdf pl_us_hazard_1hz_2pc50.ps ../figs/pl_us_hazard_1hz_2pc50.pdf
mv pl_us_hazard_1hz_2pc50.ps ../figs/
ps2pdf pl_us_hazard_5hz_2pc50.ps ../figs/pl_us_hazard_5hz_2pc50.pdf
mv pl_us_hazard_5hz_2pc50.ps ../figs/

# send mail to user
time2=`date`
mail -s "`basename $0` completed" $email_address << ENDMAIL
$0 script completed

Start time: $time1
End time: $time2

All hazard calculations written to the HERE

ENDMAIL

