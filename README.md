# gemini-load-nf
<hr/>
'''
// ******************** Start Params *********************
params.rootDir = "/home/dev/Code"
params.mapping = "${params.rootDir}/data/mapping.csv"
params.vcfs = "${params.rootDir}/data/vcfs"
params.refs = "${params.rootDir}/data/references"

params.dbs = "${params.rootDir}/data/dbs"
params.annos = "${params.rootDir}/data/annos"

params.annotateMemory = '4.2 GB'
params.loadCpus = 1
// ******************** End Params *********************
'''