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
    val  db

    output:
    path "${db}.*" 
    val "bowtie_db"

    script:
    """
    bowtie2-build $fasta $db
    """

}

process bowtie2Mapping {
    container "quay.io/biocontainers/bowtie2:2.5.4--he96a11b_5"

input:
    path fastaSubset
    val  dbname

    output:
    path "probes.sam"  
        
    script:
    """ 
    bowtie2  -x $dbname -f -U $fastaSubset -S probes.sam
    """ 
}



workflow {
    if(params.wantSplicedAlignments) {

        gmapDb = gmapBuild(params.genomeFasta)

        iit = makeKnownSpliceSiteFile(params.gtfFile)

        resultSubset = gsnapMapping(fastaSubset, iit, gmapDb)

        sam2bam(resultSubset.collectFile(name: "merged.sam"))

    }
    else {

        bowtie2Db = bowtie2Index(params.genomeFasta,params.db)

	bowtieSubset = bowtie2Mapping(fastaSubset, params.db)

	sam2bam(bowtieSubset.collectFile(name: "merged.sam"))


        //$cmd = "($bowtie2 --rg-id EuP --rg 'SM:TU114' --rg 'PL:Illumina' ";
        //if ($extraBowtieParams){$cmd = $cmd.$extraBowtieParams;}
        //$cmd = $cmd." -x $bowtieIndex ".(-e "$mateB" ? "-1 $mateA -2 $mateB " : "-U $mateA ")."-S $workingDir/$tmpOut.sam) >& $workingDir/bowtie.log";
        //sam2bam(resultSubset.collectFile(name: "merged.sam"))
    }
}

