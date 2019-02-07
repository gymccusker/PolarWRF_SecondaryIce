#!/bin/bash
#
# Automatically setup Polar WRF files
#
# Run this script from the directory 
# where you have uncompressed the WRFV3 and 
# PWRF tar files.  Change the version of WRF
# below where appropriate (search 'PWRF3.')
#  
############################################

PWRF_DIR='PWRF3.6.1'
WRF_DIR='WRFV3'
THIS_DIR=`pwd`

# Copy PWRF dir & update ${PWRF_DIR}
cp -R ${PWRF_DIR} ${WRF_DIR}/
PWRF_DIR=${WRF_DIR}/${PWRF_DIR}  

# list all Polar WRF files
PWRF_LIST=`find $PWRF_DIR/ -type f`

########################
# Loop through all files
########################

for i in $PWRF_LIST
do

	SUB_DIR=`dirname $i`
	PWRF_FILE=$SUB_DIR/$i
	PWRF_FILE_RENAMED=$SUB_DIR/`basename $i .PWRF3.6.1`
	WRF_FILE=`echo $PWRF_FILE_RENAMED | sed 's/\/PWRF3.6.1/\/./g'`

	if [ -f ${WRF_FILE} ] 
	then
		mv ${WRF_FILE} ${WRF_FILE}-unpolar
		ln -s $THIS_DIR/${PWRF_FILE_RENAMED}.PWRF3.6.1 $THIS_DIR/${WRF_FILE}
		echo `basename ${WRF_FILE}` ' : renamed old and made link to polar WRF file'
	fi

done

exit
