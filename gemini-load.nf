#!/usr/bin/env nextflow

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

def toPrefixTuple = { file -> tuple(file.name.take(file.name.lastIndexOf('.')), file.toRealPath()) }
def flatGroupTuple = { it -> tuple(it[0], *it[1].sort()) }
def fromCSVToMetas = { path -> 
  file(path).text.split("\\r?\\n").collect { 
    def vals = it.split(",")
    return [ "vcf": vals[0], "ref": vals[1], "refDB": vals[2] ]
  }
}
def resolveFile = { meta -> tuple(meta, file("${params.vcfs}/${meta.vcf}.vcf"), file("${params.refs}/${meta.ref}")}

vcfMetas = Channel.from(*fromCSVToMetas(params.mapping)).map(resolveFile)

process decomposeNormalize {

    input:
    set meta, 'file.vcf', 'ref.fasta' from vcfMetas

    output:
    set meta, "decomposed.normalized.vcf" into decompoedNormalizedVCFs

    shell: 
    '''
    zless file.vcf |
        sed 's/ID=AD,Number=./ID=AD,Number=R/' | 
        vt decompose -s - |
        vt normalize -r ref.fasta - > decomposed.normalized.vcf
    '''
}

process annotation {

    memory params.annotateMemory

    input:
    set meta, "decomposed.normalized.vcf" from decompoedNormalizedVCFs

    output:
    set meta, "annotated.vcf.gz" into annotatedVCFs
    set "snpEff_genes.txt", "snpEff_summary.html" into snpEffLogs

    shell:
    '''
    java -Xmx4G -jar $SNPEFF_JAR -v !{meta.refDB} decomposed.normalized.vcf | bgzip -c > annotated.vcf.gz
    tabix -p vcf annotated.vcf.gz
    '''
}

process geminiLoad {

    cpus params.loadCpus

    publishDir params.dbs, mode: 'copy'

    containerOptions = "-B ${params.annos}:/gemini_data"

    input: 
    set meta, "annotated.vcf.gz" from annotatedVCFs

    output: 
    set meta, "${meta.vcf}.db" into dbs

    shell:
    '''
    gemini load --cores !{params.loadCpus} -t snpEff -v annotated.vcf.gz !{meta.vcf}.db
    '''
}
    
dbs.subscribe {
    println it
}
