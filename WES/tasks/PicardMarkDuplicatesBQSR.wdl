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

    # call PicardMD {
    #     input:
    #         sample = sample,
    #         bamFile = bamFile,
    # }

    call GATK4 {
        input:
            sample = sample,
            bamFile = bamFile,
            referenceFastaFiles = referenceFastaFiles,
            thousG = thousG,
            knownIndels = knownIndels,
            dbsnp = dbsnp
    }

    output {
        File finalBam = GATK4.finalBam
    }
}

# task PicardMD {
#     input {
#         String sample
#         File bamFile
#     }

#     command <<<
#         java -Xmx1g -jar /usr/gitc/picard.jar MarkDuplicates I=~{bamFile} O=~{sample}.md.bam ASSUME_SORT_ORDER=coordinate METRICS_FILE=~{sample}.md.txt QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT
#     >>>

#     output {
#         File finalBam = "~{sample}.md.bam"
#         File metricFile = "~{sample}.md.txt"
#     }


#     runtime {
#         docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
#         disks: "local-disk 100 SSD"
#         memory: "16G"
#         cpu: 2
#     }
# }

task GATK4 {
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
        /usr/gitc/gatk4/gatk-launch BaseRecalibrator -R ~{referenceFasta} -I ~{bamFile} -O ~{sample}.bqsr.table --knownSites ~{thousG} --knownSites ~{knownIndels} --knownSites ~{dbsnp}

        /usr/gitc/gatk4/gatk-launch ApplyBQSR -R ~{referenceFasta} -I ~{bamFile} -O ~{sample}.FINAL.bam -bqsr ~{sample}.bqsr.table --static_quantized_quals 10 --static_quantized_quals 20 --static_quantized_quals 30
    >>>

    output {
        File finalBam = "~{sample}.FINAL.bam"
    }

    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        disks: "local-disk 100 SSD"
        memory: "16G"
        cpu: 2
    }
} 