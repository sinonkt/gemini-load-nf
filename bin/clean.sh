ROOT_DIR=/home/dev/Code
mv ${ROOT_DIR}/works/singularity $ROOT_DIR
rm -rf ${ROOT_DIR}/.nextflow/*
rm -rf ${ROOT_DIR}/logs/*
rm -rf ${ROOT_DIR}/works/*
mv ${ROOT_DIR}/singularity ${ROOT_DIR}/works