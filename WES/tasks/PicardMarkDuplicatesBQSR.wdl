version 1.0

workflow GatkCommands {

    input {
        String sample
        File bamFile
        Array[File] referenceFastaFiles
        Array[File] referenceGATK4Files
    }

    call PicardMD {
        input:
            sample = sample,
            bamFile = bamFile,
            referenceFastaFiles = referenceFastaFiles,
            referenceGATK4Files = referenceGATK4Files        
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
        Array[File] referenceGATK4Files
    }

    File referenceFasta = referenceFastaFiles[0]
    File knownIndels = referenceGATK4Files[0]
    File thousG = referenceGATK4Files[2]
    File dbsnp = referenceGATK4Files[4]
    
    command <<<
        java -Xmx1g -jar /usr/gitc/picard.jar MarkDuplicates I=~{bamFile} O=~{sample}.md.bam ASSUME_SORT_ORDER=coordinate METRICS_FILE=~{sample}.md.txt QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT
        
        /usr/gitc/gatk4/gatk-launch BaseRecalibrator -R ~{referenceFasta} -I ~{sample}.md.bam -O ~{sample}.bqsr.table -knownSites ~{thousG} -knownSites ~{knownIndels} -knownSites ~{dbsnp}

        /usr/gitc/gatk4/gatk-launch ApplyBQSR -R ~{referenceFasta} -I ~{sample}.md.bam -O ~{sample}.FINAL.bam -bqsr ~{sample}.bqsr.table --static_quantized_quals 10 --static_quantized_quals 20 --static_quantized_quals 30  
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