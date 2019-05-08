version 1.0

task mergeBams {
    input {
        String sample
        Array[File] bamFiles
    }

    command <<<
        samtools merge ~{sample}.bam ~{sep=" " bamFiles}
    >>>

    output {
        File mergedBam = "~{sample}.bam"
    }

    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        disks: "local-disk 100 SSD"
        memory: "8G"
        cpu: 2
    }
}
