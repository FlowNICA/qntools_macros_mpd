#!/bin/bash

#
#$ -wd /scratch2/$USER/TMP
#$ -cwd
#$ -l h=!(ncx152.jinr.ru|ncx205.jinr.ru|ncx123.jinr.ru|ncx111.jinr.ru|ncx113.jinr.ru|ncx149.jinr.ru|ncx223.jinr.ru)
#$ -N run_correlate
#$ -q all.q
#$ -l h_rt=2:30:00
#$ -l s_rt=2:30:00
#$ -t 1-10
#
#$ -o /scratch2/$USER/TMP
#$ -e /scratch2/$USER/TMP
#

source /cvmfs/nica.jinr.ru/sw/os/login.sh
module add mpddev
source /scratch2/parfenov/Soft/mpdroot/install/config/env.sh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scratch2/parfenov/Soft/QnTools/install/lib
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:/scratch2/parfenov/Soft/QnTools/install/include/QnTools

export JOB_ID=${JOB_ID}
export TASK_ID=${SGE_TASK_ID}

export FILELIST=/scratch2/parfenov/Soft/qntools_macros_mpd/macros/qn_urqmd_bibi_2.5gev_mpdfxt.list

export SHORTNAME1=`basename $FILELIST`
export SHORTNAME11=${SHORTNAME1%.list}
export LABEL=${SHORTNAME11}

export INFILE=`sed "${TASK_ID}q;d" ${FILELIST}`

if [[ -f "$INFILE" ]]; then
export DATE=${JOB_ID} # or `date '+%Y%m%d_%H%M%S'`

export MAIN_DIR=/scratch2/parfenov/Soft/qntools_macros_mpd
export OUT_DIR=${MAIN_DIR}/OUT/corr_${LABEL}/${DATE}
export OUT_FILE_DIR=${OUT_DIR}/files
export OUT_LOG_DIR=${OUT_DIR}/log
export OUT_FILE=${OUT_FILE_DIR}/corr_${LABEL}_${JOB_ID}_${TASK_ID}.root
export OUT_LOG=${OUT_LOG_DIR}/qvect_${LABEL}_${JOB_ID}_${TASK_ID}.log

export TMP=${MAIN_DIR}/TMP
export TMP_DIR=${TMP}/TMP_${JOB_ID}_${TASK_ID}

mkdir -p $OUT_FILE_DIR
mkdir -p $OUT_LOG_DIR
#mkdir -p $TMP_DIR
touch $OUT_LOG

export MACRO_EXE=${MAIN_DIR}/correlate.C

# Main process
echo "Job Id:  ${JOB_ID}" &>> $OUT_LOG
echo "Task Id: ${TASK_ID}" &>> $OUT_LOG
echo "INFILE:  ${INFILE}" &>> $OUT_LOG
echo "OUT FILE:${OUT_FILE}" &>> $OUT_LOG
echo "CONFIG:  ${MACRO_EXE}" &>> $OUT_LOG
echo "Config ROOT-macro contents:" &>> $OUT_LOG
echo "------------------------------------------------------------------------------" &>> $OUT_LOG
echo "" &>> $OUT_LOG
cat $MACRO_EXE &>> $OUT_LOG
echo "" &>> $OUT_LOG
echo "------------------------------------------------------------------------------" &>> $OUT_LOG
echo "" &>> $OUT_LOG

#cd $TMP_DIR

root -l -b -q ${MACRO_EXE}'("'${INFILE}'","'${OUT_FILE}'")' &>> $OUT_LOG

echo "Job is finished!" &>> $OUT_LOG

fi
