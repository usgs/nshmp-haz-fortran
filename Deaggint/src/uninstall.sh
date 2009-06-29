#!/bin/bash

for file in `ls`; do
	if [ -h $file ]; then
		unlink $file
	fi
done

exit 0;
