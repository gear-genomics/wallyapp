#!/bin/bash

if [ $# -lt 4 ]
then
    echo "**********************************************************************"
    echo "Wally Web Application: This program comes with ABSOLUTELY NO WARRANTY."
    echo "Version: 0.1.2"
    echo "**********************************************************************"
    echo ""
    echo "Usage: $0 <outdir> <dataset> <chr8:12-14> <NA12878> ..."
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# Terminate on error
set -e

OUTDIR=$1
DATASET=$2

## Check reference
REF=`echo ${DATASET} | sed 's/^.*_//'`
if [ ! -f ${BASEDIR}/genome/${REF}.fa.gz ]
then
    >&2 echo "ERROR: Reference genome is missing!" ${REF}.fa.gz
    exit 1;
fi

## Check output directory
if [ ! -d ${OUTDIR} ]
then
    >&2 echo "ERROR: Output directory does not exist: " ${OUTDIR}
    exit 1;
fi

## Check region
cd ${OUTDIR}
REGION=$3
if [ `echo ${REGION} | grep -c -P "^[A-Za-z0-9]*:[0-9]*-[0-9]*:[A-Za-z0-9\-]*$"` -ne 1 ]
then
    >&2 echo "ERROR: Incorrect region format: " ${REGION}
    exit 1;
fi

## Slice BAMs
shift 3
for SAMPLE in $@
do
    if [ `echo ${SAMPLE} | awk '{print length($1);}'` -eq 7 ]
    then
	REGSHORT=`echo ${REGION} | sed 's/:[^:]*$//'`
	URL=`grep -w "${SAMPLE}" ${BASEDIR}/datasets/${DATASET}.tsv`
	if [ $? -ne 0 ]
	then
	    >&2 echo "ERROR: Sample is not available in the data set: " ${SAMPLE}
	    exit 1;
	fi
	samtools view --reference ${BASEDIR}/genome/${REF}.fa.gz -b ${URL} ${REGSHORT} > ${SAMPLE}.bam
	if [ $? -ne 0 ]
	then
	    >&2 echo "ERROR: Samtools failed!"
	    exit 1;
	fi
	samtools index ${SAMPLE}.bam
	if [ $? -ne 0 ]
	then
	    >&2 echo "ERROR: Samtools indexing failed!"
	    exit 1;
	fi
    fi
done

## Plot
NUMS=`ls *.bam 2> /dev/null | wc -l`
if [ ${NUMS} -gt 0 ]
then
    HEIGHT=`echo "${NUMS} * 360" | bc -l`
    wally region -y ${HEIGHT} -pcu -r ${REGION} -g ${BASEDIR}/genome/${REF}.fa.gz *.bam
    if [ $? -ne 0 ]
    then
	>&2 echo "ERROR: Wally command-line application failed!"
	exit 1;
    fi
    rm -f *.bam *.bam.bai *.cram *.cram.crai
else
    if [ $? -ne 0 ]
    then
	>&2 echo "ERROR: No input alignment files present!"
	exit 1;
    fi
fi
