#!/bin/csh
# CSH script to run 2008-update PSH deaggregations interactively, for 48 states USA
# 
# Differs from run_2008 by including option to deagg by GMPE as well as mean haz
# $1 = pid
# $2 = output directory name
# $3 = site name no blanks, 1 <=name length <= 16 chars
# If script is passed null site_name, args 5,6,... shift left & Script becomes confused.
# $4 = site latitude, decimal degrees
# $5 = site longitude, decimal degrees
# $6 = Probability 1, 2, 5, 10, 20, 50  (new 2008)
# $7 = Number of years: 30, 50, 75, 100, 200 (new 2008)
#  
# $8 = SA frequency: 1=pga, 2=10hz, 3=5.0 hz, 4=3.33 hz, 5=2 hz, 6=1 hz 7=2.0 sec,8= 3sec. 9=4 sec, 10 = 5 sec
#     3 to 5-s are only available west of the Rockies. other frequencies ok everywhere.
# $9 = Vs30 (which implies NEHRP Site Class: A BC C CD or D (5 possibilities) WUS;
#    A or BC for CEUS; class A is 2000 m/s; BC is 760 m/s
# $10 = 0 or 1. O for mean hazard only. 1 to also output deaggs for each GMPE
#		Arg 10 is what is new in 2009 version of this cshell script.
# $11 = geographic deaggregations: 1 or 2 =run with coarse distance bin 
#				  3 or 4 = run with fine distance bin. 1&3 coarse angle bin.
#				   0 no geographic
# $12 = indicator for generating stoch. 'grams: 1=use modal R,M and 1corner SS; 
#			2=>use mean R,M, and 1-corner Source Spectrum;
#			3 = use mode and use 2-corner source spectrum. 4=use mean and 2-corSS
#			0=don't generate Stochastic seismograms.
# Other notes:
# File names use period or SA frequency. File names now unique by including PID
# A key requirement is no "standard" or "error" output. 
#  When script attempts file deletion (rm...), if not present a "no match" message occurs. This
# message is redirected to /dev/null below. This has been a big stumbling block.
# 2002 revision: now geographic deaggs can be run for all of the SA periods.
# Two directories are important: input and gmt scripts in $2 (some subdirs as well).
# Also a LOG subdirectory under $2 must exist for storing error & other diagnostics 
# Most temporary output is written to $2. $2 is periodically cleaned up.
# Otherwise, aborted runs may result in pileup of junk in $2.
set per=(0.00 0.10 0.20 0.30 0.50 1.00 2.00 3.00 4.00 5.00 )
# changed 2008 so that each extension has 3 chars. this will probably mess up
# some of the earlier stuff
cd $2
set fr= (pga 10h 5hz 3hz 2hz 1hz 2sc 3sc 4sc 5sc)
set frq= (pga 10hz 5hz 3hz 2hz 1hz 2sc 3sc 4sc 5sc)
set gmt = '/usr/local/gmt/bin/'
set pid = $1
set AZ = $11
if  ($11 > 0 ) then
set geo = 1
echo "Geog deagg was requested" >>& $2LOG/$pid.log
else
set geo = 0
echo "Geog deagg was not requested" >>& $2LOG/$pid.log
endif
# SG below is set to $12, becomes "true" when the seismogram option is exercised.
set SG = $11
if ( $8 < 1 || $8 > 10 ) then
echo "Invalid frequency flag in run_2008" >>& $2LOG/$pid.log
exit (19)
endif
if ( $8 > 7 && $5:r > -110 ) then
echo "3s and longer T SA are not available east of 110 degrees" >>& $2LOG/$pid.log
exit 19
endif
if ( $4:r < 31 && $5:r < -110 && $5:r > -154 ) then
# Pacific Ocean SW of Calif: no calcs
exit (18)
else if ( $4:r < 30 && $5:r < -105 && $5:r > -154 ) then
# Mexico or Pacific Ocean
exit (17)
else if ( $4:r < 28 && $5:r < -100 && $5:r > -154 ) then
# Mexico
exit (17)
else if ( $4:r < 28 && $5:r < -96 && $5:r > -83) then
# Gulf of Mexico
exit (16)
else if ( $4:r < 25  && $5:r > -154 ) then
exit ( 17)
else if ( $4:r == 25 && $5:r < -75 && $5:r > -125) then
exit ( 16)
else if ($4:r < 32 && $5:r > -80) then
exit (15)
# Atlantic Ocean, vicinity of Bahamas
else if ($4:r < 40 && $5:r > -74) then
exit (15)
# Atlantic Ocean, east of Carolinas etc
else if ($5:r > -66) then
exit (15)
#atlantic ocean
else if ($4:r > 45 && $5:r > -82 && $5:r < -73) then
exit (14)
# Canada (this is a crude block, lots of major cities remain inbounds)
else if ($4:r > 49 && $5:r > -131) then
exit (13)
# Canada or Atlantic Ocean
else if ($4:r > 50 && $5:r > -131) then
# Canada
exit (13)
else if ($4:r < 56 && $5:r < -139 && $5:r > -154) then
# Pacific Ocean
exit (12)
else if ( $5:r < -171) then
# Pacific Ocean or Western aleutians.. could calculate but dont
exit (12)
else if ($5:r < -122 && $5:r > -154 && $4:r < 36) then
# Pacific Ocean west of California and east of Hawaii
exit (18)
else if ($5:r < -123 && $5:r > -154 && $4:r < 37) then
# Pacific Ocean west of California and east of Hawaii
exit (18)
endif
# Below part of script is for conterminous US. No longer have 4p or 3p. just 1p.
# There is a new Alaska model but this has not been put together yet (8/2008).
set ann =  `  echo $6 $7 | awk '{print -log(1.-0.01*$1)/$2}' `
# echo $ann
$2bin/combine.2009 $3 $4 $5 $ann $8 $1  $2 $2 $9 $10 $SG $AZ  >>& $2LOG/$pid.log
set BC = `$2bin/vancouver $4 $5 `
chmod +x $2$pid.1p
#cp $2$pid.1p $2$pid.1p.sav
# Expect a log file to exist under $2 (this is a 2008 assumption)
$2$pid.1p $2 $geo $2 >>& $2LOG/$pid.log
# when it is working /bin/rm $2$pid.1p
# Currently $geo just a placeholder for e.
echo "*** Deaggregation of Seismic Hazard at One Period of Spectral Accel. ***" > $2$pid.txt
echo "*** Data from U.S.G.S. National Seismic Hazards Mapping Project, 2008 version ***" >> $2$pid.txt
#echo "*** Spectral periods 0.1 s, 0.3 s, 0.5 s, and 2.0 s are coming soon! *****" >> $2$pid.txt
# file names in 2008 are PID // 3-letter source-model abbrev. Remove when process is well-tested
/bin/rm $2$pid??? >& /dev/null
/bin/rm $2$pid*.FDEAG >& /dev/null
# the above removes unwanted temp files (there may be some by mistake)

#echo the .1p file finished running
cat $2$pid.$fr[$8] >> $2$pid.txt
if ($BC == 1) then
echo "******************** Vancouver Island region ************************************" >> $2$pid.txt
else if ($BC == 2) then
echo "******************** Pacific Northwest region ************************************" >> $2$pid.txt
else if ($BC == o) then
echo "******************** Southern Oregon site ************************************" >> $2$pid.txt
else if ($BC == 3) then
echo "******************** Southern California ****************************************" >> $2$pid.txt
else if ($BC == 4) then
echo "******************** Central  California ****************************************" >> $2$pid.txt
else if ($BC == 5) then
echo "******************** Northern California ****************************************" >> $2$pid.txt
else if ($BC == N) then
echo "******************** Northern Nevada Site *************************************" >> $2$pid.txt
else if ($BC == 6) then
echo "******************** Intermountain Seismic Belt***********************************" >> $2$pid.txt
else if ($BC == 7) then
echo "******************** Basin and Range site  *************************************" >> $2$pid.txt
else if ($BC == 8 ) then
echo "********************Central or Eastern U.S. Site ********************************" >> $2$pid.txt
else if ($BC == 9) then
echo "******************** Southern Canada Site  **************************************" >> $2$pid.txt
else if ($BC == C) then
echo "******************** Colorado Site *********************************************" >> $2$pid.txt
else if ($BC == T) then
echo "******************** Texas Site *********************************************" >> $2$pid.txt
else if ($BC == O) then
echo "******************** Oklahoma Site *********************************************" >> $2$pid.txt
else if ($BC == S ) then
echo "******************** South Carolina & Vicinity****************************************" >> $2$pid.txt
else if ($BC == M) then
echo "******************** Mexico Site *********************************************" >> $2$pid.txt
else
echo "********************************************************************************" >> $2$pid.txt
endif
/bin/rm $2$pid.$fr[$8] >& /dev/null
# The file $pid.txt could be moved to ftp area. It contains src deagg for one SA frequency
chmod +x $2gmt.$per[$8]$pid
$2gmt.$per[$8]$pid $gmt  $2 >>& $2LOG/$pid.log
if (-e $2$pid$per[$8].ps ) mv $2$pid$per[$8].ps $2$pid$frq[$8].ps >& /dev/null
# Output postscript is $2$pid$frq[$8].ps for this case
if ($SG && -e $2smsim.gmt.$pid ) then
chmod +x $2smsim.gmt.$pid
$2smsim.gmt.$pid $gmt $2 >>& $2LOG/$pid.log
endif
if ($geo && -e $2gmtg.$per[$8].$pid ) then
#/usr/bin/cp $2/PL/not.ready.ps $2$pid$frq[$8].ps
	chmod +x $2gmtg.$per[$8].$pid
	$2gmtg.$per[$8].$pid $GMTN $2 >>& $2LOG/$pid.log
# 2002 version: wwwdefaults and others are located in the subdirectory PL
if ( -e $2$pid$per[$8].g.ps )	mv $2$pid$per[$8].g.ps $2$pid$frq[$8]g.ps
# Geographic plot is $2$pid$frq[$8]g.ps for this case
# Geographic scripts and cpt files are removed next
	/bin/rm $2$per[$8]$pid.cpt >& /dev/null
	/bin/rm $2gmtg.$per[$8].$pid >& /dev/null
# jan 2003 geog data files have suffix .geo
	/bin/rm $2$pid*.?0.geo >& /dev/null
endif
# The file $pid.txt will be moved to ftp area. It contains SA for 4 frequencies.
/bin/rm $2addhaz.$pid*  >& /dev/null
/bin/rm $2addhaz.${pid}A  >& /dev/null
/bin/rm $2gmt.*$pid >& /dev/null
/bin/rm $2gmset.*$pid >& /dev/null
if ($SG) then
	 /bin/rm $2smsim.gmt.$pid 
endif
exit (0)
