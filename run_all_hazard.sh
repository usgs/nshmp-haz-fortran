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
  echo "USAGE: $0 [path to NSHMP directory] [email address for final notification]"
  echo "Email will be sent to user when processing is completed"
  exit
fi
trunk_path=$1
email_address=$2

# 
scripts_path=$trunk_path/scripts
time1=`date`

# change to scripts directory
# check all scripts are executable
if [ -d $scripts_path ]; then
  cd $scripts_path
  chmod 750 hazrun_*.sh
  chmod 750 combine2014_all.sh
else
  echo "$scripts path does not exist"
  exit
fi
 
# create the out and logs directories for output 
if [ ! -d ../out ]; then
  mkdir ../out
  mkdir ../out/combine
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

# create all required rjb files
cd ../bin
if [ ! -f meanrjb.bin ]; then
./getmeanrjf.v2 << END
1
END
  mv rjbmean.bin.srl meanrjb.bin
  mv rjbmean.dat.srl meanrjb_srl.dat
./getmeanrjf.v2 << END
2
END
./getmeanrjf.v2 << END
3
END
fi

# cd to scripts directory for consistency with hazard input files
cd $scripts_path

# run hazard calculation scripts
for scr_nm in hazrun_casc_2014.sh hazrun_ceus_2014.sh hazrun_wus_2014.sh; do
  reg=`echo $scr_nm | cut -d_ -f2`
  logf=../logs/log_${reg}_2014.txt
  echo "Running $scr_nm - $logf"
  ./$scr_nm >& $logf
done
time5a=`date`

# send mail to user
time2=`date`
mail -s "`basename $0` completed" $email_address << ENDMAIL
$0 script completed

Start time: $time1
End time: $time2

NO Plotting in this version of the script


All hazard calculations written to files in out/ directory.

ENDMAIL
