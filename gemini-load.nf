#!/usr/bin/env nextflow

// ******************** Start Params *********************
params.rootDir = "/home/dev/Code"
params.tempDir = "${params.rootDir}/temp"
params.mapping = "${params.rootDir}/data/mapping.csv"
params.vcfs = "${params.rootDir}/data/vcfs"
params.refs = "${params.rootDir}/data/references"

params.dbs = "${params.rootDir}/data/dbs"
params.annos = "${params.rootDir}/data/annos"

params.annotateMemoryGB = 4
params.annotateCpus = 1
params.loadCpus = 1
// ******************** End Params *********************

def fromCSVToMetas = { path -> 
  file(path).text.split("\\r?\\n").collect { 
    def vals = it.split(",")
    splitted = vals[0].tokenize('.')
    vcfFileId = splitted.take(splitted.size()-2).join(".") // remove .vcf.gz
    return [ "vcf": vcfFileId, "ref": vals[1], "refDB": vals[2], "ped": vals.length == 4 ? vals[3]: null ]
  }
}

def resolveFile = { meta -> tuple(meta, 
    file("${params.vcfs}/${meta.vcf}.vcf.gz"),
    file("${params.refs}/${meta.ref}"),
    file("${params.refs}/${meta.ped}"))
}

vcfMetas = Channel.from(*fromCSVToMetas(params.mapping)).map(resolveFile)

process decomposeNormalizeAnnotate {

    memory params.annotateMemoryGB

    input:
    set meta, 'file.vcf.gz', 'ref.fasta' from vcfMetas

    output:
    set meta, "annotated.vcf.gz", "annotated.vcf.gz.tbi" into annotatedVCFs

    shell:
    '''
    bgzip --decompress --threads !{params.annotateCpus} -c file.vcf.gz |
        sed 's/ID=AD,Number=./ID=AD,Number=R/' | 
        vt decompose -s - |
        vt normalize -r ref.fasta - |
        java -Xmx!{params.annotateMemoryGB}G -jar $SNPEFF_JAR -t !{meta.refDB} |
        bgzip --threads !{params.annotateCpus} -c > annotated.vcf.gz
    tabix -p vcf annotated.vcf.gz
    '''
}

process geminiLoad {

    publishDir params.dbs, mode: 'copy'

    containerOptions = "-B ${params.annos}:/gemini_data"

    input: 
    set meta, "annotated.vcf.gz", "annotated.vcf.gz.tbi" from annotatedVCFs

    output: 
    set meta, "${meta.vcf}.db" into dbs

    shell:
    if(meta.ped == null)
        '''
        gemini load --cores !{params.loadCpus} --tempdir !{params.tempDir} -t snpEff -v annotated.vcf.gz \
            !{meta.vcf}.db
        '''
    else
        '''
        gemini load --cores !{params.loadCpus} --tempdir !{params.tempDir} -t snpEff -v annotated.vcf.gz \
            -p !{meta.pedFile} \
            !{meta.vcf}.db
        '''
}
    
dbs.subscribe {
    println it
}
