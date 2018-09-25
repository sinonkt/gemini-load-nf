ROOT_DIR=/home/dev/Code
datetime=$(date -d "today" +"%Y%m%d%H%M")
NX_CONFIG=nextflow.config
LOG_DIR=${ROOT_DIR}/logs/${datetime}
WORK_DIR=${ROOT_DIR}/works/

mkdir -p $LOG_DIR

nextflow -log $LOG_DIR/.nextflow.log \
  -C $NX_CONFIG \
  run $1 \
  -w $WORK_DIR \
  -resume \
  -with-report ${LOG_DIR}/report.html \
  -with-dag ${LOG_DIR}/flowchart.svg \
  -with-timeline ${LOG_DIR}/timeline.html \
