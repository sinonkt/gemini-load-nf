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

vcfMetas = Channel.from(*fromCSVToMetas(params.mapping))
vcfs = Channel.fromPath("${params.vcfs}/*.vcf")
references = Channel.fromPath("${params.refs}/*")


process decomposeNormalize {

    input:
    val meta from vcfMetas
    file "${meta.vcf}.vcf" from vcfs
    file "${meta.ref}" from references

    output:
    set meta, "decomposed.normalized.vcf" into decompoedNormalizedVCFs

    shell: 
    '''
    zless !{meta.vcf}.vcf |
        sed 's/ID=AD,Number=./ID=AD,Number=R/' | 
        vt decompose -s - |
        vt normalize -r !{meta.ref} - > decomposed.normalized.vcf
    '''
}

process annotation {

    memory params.annotateMemory

    input:
    set meta, "decomposed.normalized.vcf" from decompoedNormalizedVCFs

    output:
    set meta, "annotated.vcf.gz", "snpEff_genes.txt", "snpEff_summary.html"  into annotatedVCFs

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
    set meta, "annotated.vcf.gz", "snpEff_genes.txt", "snpEff_summary.html"  from annotatedVCFs

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
