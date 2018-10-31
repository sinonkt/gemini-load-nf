# gemini-load-nf
Implicit value of params which synchronize to [docker-centos7-singularity-nextflow](https://hub.docker.com/r/sinonkt/docker-centos7-singularity-nextflow/~/dockerfile/)
```groovy
params.rootDir = "/home/dev/Code"
params.mapping = "${params.rootDir}/data/mapping.csv"
params.vcfs = "${params.rootDir}/data/vcfs"
params.refs = "${params.rootDir}/data/references"

params.dbs = "${params.rootDir}/data/dbs"
params.annos = "${params.rootDir}/data/annos"

params.annotateMemoryGB = 4
params.annotateCpus = 1
params.loadCpus = 1
```
