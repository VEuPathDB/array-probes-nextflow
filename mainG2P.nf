#!/usr/bin/env nextflow
nextflow.enable.dsl=2

if(!params.fastaSubsetSize) {
  throw new Exception("Missing params.fastaSubsetSize")
}

if(params.probesFastaFile) {
  fastaSubset = Channel.fromPath( params.probesFastaFile )
           .splitFasta( by:params.fastaSubsetSize, file:true  )
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

    input:
      path probes.bam
      path genome.gff

    output:
      path "probes.bed"

    script:
    """
    bedtools intersect -wa -wb -a $probes.bam -b $genome.gff -bed > probes.bed
    """
}


process bed2Tab {
  
    input:
      path probes.bed

    output:
      path "gene2probe.list"


    script:
    """
      grep "ID=exon_" $probes.bed > exons.out
      cut -f4 exons.out > probes.txt
      cut -d';' -f3 exons.out > geneID.txt
      sed 's/gene_id=//g' geneID.txt > geneIDs.txt
      paste geneIDs.txt probes.txt > gene2probe.list
    """
}

process hashGene2probe {

    input:
      path gene2probe
      path gene2probe.py

    output:
      path "gene2probe.tsv"

    script:
    """
    python $gene2probe.py ($gene2probe) > gene2probe.tsv


workflow {
    if(params.wantSplicedAlignments) {

        gmapDb = gmapBuild(params.genomeFasta)

        iit = makeKnownSpliceSiteFile(params.gtfFile)

        resultSubset = gsnapMapping(fastaSubset, iit, gmapDb)

        sam2bam(resultSubset.collectFile(name: "merged.sam"))

        // IF Stranded (params.makeCdf == true) NOTE:  may need to play around with the htseq commands
        //htseq-count -a 0 --format=bam --order=name --stranded=yes --type=exon --idattr=gene_id --mode=union BAM_FILE params.gtfFile
        // ELSE
        // htseq-count -a 0 --format=bam --order=name --stranded=no --type=exon --idattr=gene_id --mode=union BAM_FILE .bam params.gtfFile

        // TODO: geneToProbeMapping (output geneToProbeMapping file)

//        if(params.makeCdfFile) {
            //$MIN_PROBES=3
            //my $cmd1 = "makePbaseTbase.pl probes.fsa $workflowDataDir/pbase-tbase.out";
            //my $cmd2 = "makeCdfHeader.pl  --outPutFile $outputCdfFile --gene2probes $gene2probesFile --name $name --rows params.arrayRows --cols  params.arrayColumns --minProbes $MIN_PROBES";
            //my $cmd3 = "create_cdf.pl $workflowDataDir/$outputCdfFile $workflowDataDir/$gene2probesInputFile $workflowDataDir/pbase-tbase.out $MIN_PROBES";
//        }
//        if(params.makeNdfFile) {
            // my $cmd = "recreate_ndf.pl --original_ndf_file $workflowDataDir/$ndfFile --gene_to_oligo_file $workflowDataDir/$gene2probesInputFile --output_file $workflowDataDir/$outputFile";
//        }

    }
    else {

        indexName = "index"

        bowtie2Db = bowtie2Index(params.genomeFasta, indexName)

        bowtieSubset = bowtie2Mapping(fastaSubset, bowtie2Db, indexName)

        sam2bam(bowtieSubset.collectFile(name: "merged.sam"))
    }
}

