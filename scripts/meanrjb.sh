#!/bin/bash

# This script will (re)create all required rjb files
# on the assumption that bin/getmeanrjf.v2 exists.
# Can be called from scripts directory of project root (make).

PATHDIR=`pwd`
DIR=`basename "$PATHDIR"`
if [[ $DIR == 'scripts' ]]; then
  echo 'yo'
  cd ..
fi
cd bin

./getmeanrjf.v2 << END
1
END
mv rjbmean.bin.srl meanrjb.bin
rm rjbmean.dat.srl

./getmeanrjf.v2 << END
3
END
rm rjbmean.dat.Aeast

./getmeanrjf.v2 << END
4
END
rm rjbmean.dat.geoma
