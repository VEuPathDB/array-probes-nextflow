#!/usr/bin/env nextflow
nextflow.enable.dsl=2


// we are not spllitting probes here... to avoid merging complications due to sam file headers
if(params.probesFastaFile) {
  fastaSubset = Channel.fromPath( params.probesFastaFile )
}
else {
  throw new Exception("Missing params.probesFastaFile")
}


process gmapBuild {

    container "quay.io/biocontainers/gmap:2024.11.20--pl5321hb1d24b7_1"

    input:
    path fasta

    output:
    path "gmap_db"
    val  "gmap_db"

    script:
    """
    gmap_build -d gmap_db -D . $fasta
    """
}


process makeKnownSpliceSiteFile {
    container "quay.io/biocontainers/gmap:2024.11.20--pl5321hb1d24b7_1"

    input:
    path gtf

    output:
    path "splice_sites.iit"

    script:
    """
    gtf_splicesites $gtf | iit_store -o splice_sites.iit
    """
}


process gsnapMapping {
    container "quay.io/biocontainers/gmap:2024.11.20--pl5321hb1d24b7_1"

    input:
    path fastaSubset
    path spliceSiteFile
    path gmapDatabase
    val databaseName

    output:
    path "subset.sam"

    script:
    """
    gsnap --force-xs-dir --quiet-if-excessive -N 1 -s $spliceSiteFile  -A sam -n ${task.ext.nPaths} -D $gmapDatabase -d $databaseName  $fastaSubset >subset.sam
    """
}



process sam2bam {
    publishDir params.outputDir, mode: 'copy'

    container "quay.io/biocontainers/samtools:1.20--h50ea8bc_0"

    input:
    path sam

    output:
    path params.outputFileName
    path params.outputFileName + ".bai"

    script:
    """
    samtools view -bS $sam  -o unsorted.bam
    samtools sort unsorted.bam >$params.outputFileName
    samtools index $params.outputFileName
    """

}


process bowtie2Index {
    container "quay.io/biocontainers/bowtie2:2.5.4--he96a11b_5"

    input:
    path fasta
    val  indexBaseName

    output:
    path "${indexBaseName}.*"


    script:
    """
    bowtie2-build $fasta $indexBaseName
    """

}

process bowtie2Mapping {
    container "quay.io/biocontainers/bowtie2:2.5.4--he96a11b_5"

input:
    path fastaSubset
    path indexFiles
    val  dbname

    output:
    path "probes.sam"  
        
    script:
    """ 
    bowtie2  -x $dbname -f -U $fastaSubset -S probes.sam
    """ 
}

process probeGeneIntersect {
    container "quay.io/biocontainers/bedtools:2.27.1--h077b44d_9"

    input:
    path probeAlignmentsBam
    path probeAlignmentsBai
    path genomeGff

    output:
    path "probes.bed"

    script:
    """
    bedtools intersect -wa -wb -a $probeAlignmentsBam -b $genomeGff -bed > probes.bed
    """
}

process bedToGene2Probe {
    container = 'bioperl/bioperl:stable'
    publishDir params.outputDir, mode: 'copy'

    input:
    path bed

    output:
    path params.outputMappingFileName

    script:
    """
    bedToGene2Probe.pl --gtfFeatureType ${task.ext.gtfFeatureType} --bed $bed --geneGtfTag ${task.ext.geneGtfTag} --outputFile ${params.outputMappingFileName}
    """
}


process cdfFromGene2Probe {
    container = 'bioperl/bioperl:stable'
    publishDir params.outputDir, mode: 'copy'

    input:
    path gene2probes
    path probesFasta
    val vendorPath

    output:
    path vendorPath.name

    script:
    def vendorFileName = vendorPath.name
    def vendorBaseName = vendorPath.baseName
    """
    makePbaseTbase.pl $probesFasta >pbase-tbase.out
    makeCdfHeader.pl --outPutFile $vendorFileName --gene2probes $gene2probes --name $vendorBaseName --rows ${params.arrayRows} --cols  ${params.arrayColumns} --minProbes ${task.ext.minProbes}
    create_cdf.pl $vendorFileName $gene2probes pbase-tbase.out ${task.ext.minProbes}
    """
}

process ndfFromGene2Probe {
    container = 'bioperl/bioperl:stable'
    publishDir params.outputDir, mode: 'copy'

    input:
    path gene2probes
    val vendorPath
    path stagedVendorPath, stageAs: 'input_files/vendor.ndf'

    output:
    path vendorPath.name

    script:
    def vendorFileName = vendorPath.name
    """
    recreate_ndf.pl --original_ndf_file $stagedVendorPath --gene_to_oligo_file $gene2probes --output_file $vendorFileName
    """
}

workflow {

    def vendorMappingFile = file(params.vendorMappingFile)

    if(params.wantSplicedAlignments) {

        gmapDb = gmapBuild(params.genomeFastaFile)

        iit = makeKnownSpliceSiteFile(params.gtfFile)

        resultSubset = gsnapMapping(fastaSubset, iit, gmapDb)

        sam2bam(resultSubset.collectFile(name: "merged.sam"))

        probeGeneIntersect(sam2bam.out, params.gtfFile)

        bedToGene2Probe(probeGeneIntersect.out)

        if(params.makeCdfFile) {
            cdfFromGene2Probe(bedToGene2Probe.out, params.probesFastaFile, vendorMappingFile)
        }

        if(params.makeNdfFile) {
            ndfFromGene2Probe(bedToGene2Probe.out, vendorMappingFile, params.vendorMappingFile)
        }
    }
    else {

        indexName = "index"

        bowtie2Db = bowtie2Index(params.genomeFastaFile, indexName)

        bowtieSubset = bowtie2Mapping(fastaSubset, bowtie2Db, indexName)

        sam2bam(bowtieSubset.collectFile(name: "merged.sam"))
    }
}

