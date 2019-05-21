version 1.0

workflow GatkCommands {

    input {
        String sample
        File bamFile
        Array[File] referenceFastaFiles
        File thousG
        File knownIndels
        File dbsnp
    }

    call PicardMD {
        input:
            sample = sample,
            bamFile = bamFile,
            referenceFastaFiles = referenceFastaFiles,
            thousG = thousG,
            knownIndels = knownIndels,
            dbsnp = dbsnp
    }

    output {
        File finalBam = PicardMD.finalBam
    }
}

task PicardMD {
    input {
        String sample
        File bamFile
        Array[File] referenceFastaFiles
        File thousG
        File knownIndels
        File dbsnp
    }

    File referenceFasta = referenceFastaFiles[0]

    command <<<
        java -Xmx1g -jar picard.jar MarkDuplicates I=~{bamFile} O=~{sample}.md.bam ASSUME_SORT_ORDER=coordinate METRICS_FILE=~{sample}.md.txt QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT

        /usr/gitc/gatk4/gatk-launch BaseRecalibrator -I ~{sample}.md.bam -R {referenceFasta} --known-sites ~{thousG} --known-sites ~{knownIndels} --known-sites ~{dbsnp} -O ~{sample}.bqsr.table

        /usr/gitc/gatk4/gatk-launch ApplyBQSR -R ~{referenceFasta} -I ~{sample}.md.bam -bqsr ~{sample}.bqsr.table --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 -O ~{sample}.FINAL.bam
    >>>

    output {
        File finalBam = "~{sample}.FINAL.bam"
    }


    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        disks: "local-disk 100 SSD"
        memory: "8G"
        cpu: 2
    }
}
