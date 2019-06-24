version 1.0

task IndexBams {
    input {
        String sample
        File bamFile
    }

    command <<<
        /usr/local/bin/samtools index ~{bamFile}
    >>>

    output {
        File IndexedBamFile = "~{sample}.bam.bai"
    }

    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        disks: "local-disk 100 SSD"
        memory: "16G"
        cpu: 1
    }    
}
