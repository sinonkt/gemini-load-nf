DEV_HOME=/home/dev/Code
DEV_IMAGE=sinonkt/docker-centos7-singularity-nextflow:singularity-2.6.0

docker run -it -u dev --volume $(pwd):$DEV_HOME --privileged $DEV_IMAGE /bin/bash
