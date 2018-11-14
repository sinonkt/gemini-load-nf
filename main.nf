#!/usr/bin/env nextflow
// ******************** Start Params *********************
params.inputCSV = 'input.csv' // vcf, vcf.gz, bed_bim_fam, already_by_chromosome (every format has chromosome awareness)
params.publishDir = "publish"
// ******************** End Params *********************
// ****************** Start Helpers *****************
def VALID_CHROMOSOMES=[ *(1..22), 'X', 'Y', 'MT']
def AVAILABLE_EXTENSIONS=['vcf.gz', 'vcf', 'bed', 'bim', 'fam']
def CHECK_CHROMOSOME_PATTERN = ~/.*chr(X|Y|MT|\d*).*/
ChannelByFileType = [ 
  'bed-bim-fam': [ channelIdx: 0, extPattern: ~/\.(bed|bim|fam)$/ ],
  'vcf': [ channelIdx: 1, extPattern: ~/\.vcf$/ ],
  'vcf-gz': [ channelIdx: 2, extPattern: ~/\.vcf.gz$/ ]
]
def checkDataSetFileType = { files -> 
  for (fileType in ['bed-bim-fam', 'vcf', 'vcf-gz']) {
    if (files.any { file -> file =~ ChannelByFileType[fileType].extPattern })
      return fileType
  }
  throw new Exception("Empty directory or dataset file extension mismatch."); 
}
def toDataSetInfo = {
  allFiles = file("${it.projectDir}/*").collect { it.name }
  it.datasetFileType = checkDataSetFileType(allFiles)
  it.files = allFiles.findAll({ file -> file =~ ChannelByFileType[it.datasetFileType].extPattern })
  it.fileSetId = it.files.collect({ it.replaceAll(/.bed|.bim|.fam|.vcf.gz|.vcf$/,'') }).unique()[0]
  it
}
def attachBedBimFamFile = {
  it.bedBimFamDir = it.projectDir
  it.ped = "${it.projectDir}/${it.fileSetId}.fam"
  tuple(it, *['bed', 'bim', 'fam'].collect({ ext -> file("${it.projectDir}/${it.fileSetId}.${ext}") }))
}
def attachVcfFile = { 
  it.vcfDir = it.projectDir
  it.vcf = "${it.projectDir}/${it.fileSetId}.vcf";
  tuple(it, file(it.vcf))
}
def updateProjectDir = {
  def dataset = it[0]
  def vcfgzFile = it[1]
  dataset.projectDir = vcfgzFile.parent
  dataset
}
def checkChromosomeSplittedThenAttachVcfgzFile = { dataset ->
  allFiles = file("${dataset.projectDir}/*").collect { it }
  onlyVcfgz = allFiles.findAll({ file -> file =~ ChannelByFileType['vcf-gz'].extPattern })
  dataset.allVcfgzFiles = onlyVcfgz
  allMatched = dataset.allVcfgzFiles
    .collect({ file -> file =~ CHECK_CHROMOSOME_PATTERN })
    .findAll({ matched -> matched.matches() })
  isSplitted = !allMatched.isEmpty()
  dataset.vcfgzChrPattern = isSplitted ? dataset.fileSetId.replaceAll(/\.chr.*\./, '.chr*.') : "${dataset.fileSetId}.chr*.vcf.gz"
  dataset.isSplitted = isSplitted
  dataset.chromosomes = isSplitted ?  allMatched.collect { it[0][1] } : VALID_CHROMOSOMES
  dataset.splittedChannelIdx = isSplitted ? 0 : 1  // Splitted | NoneSplitted
  dataset
}
// ******************** End Helpers ******************

BedBimFam = Channel.create()
Vcf = Channel.create()
Vcfgz = Channel.create()
Channel
    .from(file(params.inputCSV).text)
    .splitCsv(header: true)
    .map(toDataSetInfo)
    .choice (BedBimFam, Vcf, Vcfgz) { ChannelByFileType[it.datasetFileType].channelIdx }

process VCFgz {

  tag { dataset.project }

  input:
  set dataset, "${dataset.fileSetId}.bed", "${dataset.fileSetId}.bim", "${dataset.fileSetId}.fam" from BedBimFam.map(attachBedBimFamFile)

  output:
  set dataset, "${dataset.fileSetId}.vcf.gz" into convertedVcfgzFromBedBimFam

  shell:
  '''
  plink --bfile !{dataset.fileSetId} --recode vcf --out !{dataset.fileSetId}
  bgzip --threads 8 -c !{dataset.fileSetId}.vcf > !{dataset.fileSetId}.vcf.gz
  '''
}

process Bgzip {

  tag { dataset.project }

  input:
  set dataset, "${dataset.fileSetId}.vcf" from Vcf.map(attachVcfFile)

  output:
  set dataset, "${dataset.fileSetId}.vcf.gz" into convertedVcfgzFromVcf

  shell:
  '''
  bgzip --threads 8 -c !{dataset.fileSetId}.vcf > !{dataset.fileSetId}.vcf.gz
  '''
}

mixedVcfgz = Vcfgz.mix(
  convertedVcfgzFromBedBimFam.map(updateProjectDir),
  convertedVcfgzFromVcf.map(updateProjectDir)
)

SplittedVcfgz = Channel.create()
NonSplittedVcfgz = Channel.create()

mixedVcfgz
  .map(checkChromosomeSplittedThenAttachVcfgzFile)
  .choice(SplittedVcfgz, NonSplittedVcfgz) { it.splittedChannelIdx }

process splitChromosome {

  tag { dataset.project }

  input:
  set dataset, "${dataset.fileSetId}.vcf.gz" from NonSplittedVcfgz.map{ tuple(it, file(it.allVcfgzFiles[0])) }

  output:
  set dataset, "${dataset.vcfgzChrPattern}" into SplittedVcfgzFromNonSplitted

  shell:
  '''
  tabix -p vcf !{dataset.fileSetId}.vcf.gz
  for chrIdx in !{dataset.chromosomes.join(\' \')}
  do
    tabix -h !{dataset.fileSetId}.vcf.gz $chrIdx > !{dataset.fileSetId}.chr$chrIdx.vcf
    numRows=$(gawk '/^[^#]/ {print $0}' !{dataset.fileSetId}.chr$chrIdx.vcf | wc -l)
    if [ $numRows -ne 0 ]; then
      bgzip --threads 8 -c !{dataset.fileSetId}.chr$chrIdx.vcf > !{dataset.fileSetId}.chr$chrIdx.vcf.gz
    fi
  done
  '''
}

def flattenDatasetAsVCFgzChunks = {
    def out = []
    it.splittedVcfgzs.each { vcfgz ->
      def vcfgzId = vcfgz.name.replaceAll(/.vcf.gz$/,'')
      def chrIdx = (vcfgzId =~ CHECK_CHROMOSOME_PATTERN)[0][1]
      out.push([ *:it, vcfgzId: vcfgzId, chrIdx: chrIdx, vcfgzPath: vcfgz ])
    }
    out
}

SplittedVcfgz.map({ it.splittedVcfgzs = it.allVcfgzFiles; it }).mix(
  SplittedVcfgzFromNonSplitted.map { it[0].splittedVcfgzs = it.tail().flatten(); it[0] }
)
  .flatMap(flattenDatasetAsVCFgzChunks)
  .map({ tuple(it, file(it.vcfgzPath), file(it.ref)) })
  .set { DatasetChunks }

process decomposeNormalizeAnnotate {

    tag { "${chunk.project}_${chunk.chrIdx}" }

    input:
    set chunk, "file.vcf.gz", "ref.fasta" from DatasetChunks

    output:
    set chunk, "annotated.vcf.gz", "annotated.vcf.gz.tbi" into AnnotatedVCFChunks

    shell:
    '''
    bgzip --decompress --threads 8 -c file.vcf.gz |
        sed 's/ID=AD,Number=./ID=AD,Number=R/' | 
        vt decompose -s - |
        vt normalize -r ref.fasta - |
        java -Xmx8G -jar $SNPEFF_JAR -t !{chunk.refBuild} |
        bgzip --threads 8 -c > annotated.vcf.gz
    tabix -p vcf annotated.vcf.gz
    '''
}


process geminiLoad {

    tag { "${chunk.project}_${chunk.chrIdx}" }

    input: 
    set chunk, "annotated.vcf.gz", "annotated.vcf.gz.tbi" from AnnotatedVCFChunks

    output: 
    set chunk, "${chunk.vcfgzId}.db" into DBChunks

    shell:
    if(chunk.ped == null)
        '''
        gemini load --cores 32 --tempdir ./tmp -t snpEff -v annotated.vcf.gz !{chunk.vcfgzId}.db
        '''
    else
        '''
        gemini load --cores 32 --tempdir ./temp -t snpEff -v annotated.vcf.gz -p !{chunk.ped} !{chunk.vcfgzId}.db
        '''
}
    
DBChunks.subscribe {
    println "${it.project}_${it.chrIdx}"
}

// assumption each project can have multiple fileSet. but what's about ped file
// may be ^-- it's not gonna work.
// Hey you!
// findAll to filterout not interesed files.
// then match unique fileSetId  
// process splitChromosome (every things after this already by chromosome)
// process decomposeNormalizeAnnotate 
// process geminiLoad
// merge multiple way => merge by chromosome, merge by dataset                  , merge all
//                              |                   |                                |
//                              v                   v                                v
//                       all set over a chr  | all chromosome on this set | just all 
// merge db
// map && collect is the most frequenly bug happened
