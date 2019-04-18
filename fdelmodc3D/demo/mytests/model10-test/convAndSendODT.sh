#!/bin/bash

# Check inputs
if [ "$#" -ne 1 ]; then
	echo "convAndSendODT: wrong number of input parameters. Exiting."
	echo "Converts SU file to SEGY file, move it to $ODTAREA"
	echo "Usage: "
	echo "./convAndSendODT file.su"
	#return #for function
	exit #for script
fi


su2sgy.sh $1 

name="${1%%.*}" # strip name

mv $name.sgy $ODTAREA
