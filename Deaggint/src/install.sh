#!/bin/bash

while [ 1 ]; do
	read -p 'Enter location of data files: ' datadir
	if [[ ! -z $datadir && -d "${datadir}/GM" ]]; then
		break;
	fi
	echo "The directory:"
	echo "    ${datadir}"
	echo -n "does not exist or does not or does not appear to contain "
	echo "appropriate datafiles."
done

for file in `ls $datadir`; do
	if [ -h $file ]; then
		unlink $file;
	fi
	ln -sf ${datadir}/${file} ${file}
done

exit 0;
