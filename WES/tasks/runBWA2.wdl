version 1.0

workflow runBWA {
    input {
        String sample
        Int j
        File fastq1
        File fastq2
        File referenceFasta 
    }

    call BWACommand {
            input:
                j = j,
                sample = sample,
                fastq1 = fastq1,
                fastq2 = fastq2,
                referenceFasta = referenceFasta,
        }

    output {
        Array[String] Look = BWACommand.Look
    }
}

task BWACommand {
    input {
        Int j
        String sample
        File fastq1
        File fastq2
        File referenceFasta
    }

    command <<<
        id=~(zcat < ~{fastq1} | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/./g')
        pu=~(zcat < ~{fastq1} | head -n 1 | cut -f 3-4 -d":" | sed 's/@//' | sed 's/:/./g')
        sm=~(zcat < ~{fastq1} | head -n 1 | grep -Eo "[ATGCN]+$")        
        echo "/usr/gitc/bwa mem -K 100000000 -t 8 -R @RG\tID:~id\tPU:~pu.~sm\tSM:~{sample}\tLB:~id.~sm\tPL:ILLUMINA\tCN:UCLA ~{referenceFasta} ~{fastq1} ~{fastq2} > ~{sample}.~{j}.sam" >> Look
    >>>

    output {
        Array[String] Look = read_lines("Look")
    }

    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        disks: "local-disk 100 SSD"
        memory: "8G"
        cpu: 2
    }
}
