version 1.0

workflow runBWA {
    input {
        String sample
        Array[Array[String]] fileList #tsv of sample, fastq1/2, index
        Array[File] referenceFastaFiles
    }

    scatter (filePair in fileList) {
        String sample = filePair[0]
        File fastq1 = filePair[1]
        File fastq2 = filePair[2]
        Int j = filePair[3]

        call BWACommand {
            input:
                j = j,
                sample = sample,
                fastq1 = fastq1,
                fastq2 = fastq2,
                referenceFastaFiles = referenceFastaFiles,
        }

        call Samblaster {
            input:
                sample = sample,
                samFile = BWACommand.samFile
        }

        call toBam {
            input:
                j = j,
                sample = sample,
                blastsamFile = Samblaster.blastsamFile
        }
    }

    output {
        Array[File] bamFile = toBam.bamFile
    }
}

task BWACommand {
    input {
        Int j
        String sample
        File fastq1
        File fastq2
        Array[File]+ referenceFastaFiles
    }

    File referenceFasta = referenceFastaFiles[0]

    command <<<
        id=$(zcat < ~{fastq1} | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/./g')
        pu=$(zcat < ~{fastq1} | head -n 1 | cut -f 3-4 -d":" | sed 's/@//' | sed 's/:/./g')
        sm=$(zcat < ~{fastq1} | head -n 1 | grep -Eo "[ATGCN]+$")
        /usr/gitc/bwa mem -K 100000000 -t 8 -R "@RG\tID:$id\tPU:$pu.$sm\tSM:~{sample}\tLB:$id.sm\tPL:ILLUMINA\tCN:UCLA" ~{referenceFasta} ~{fastq1} ~{fastq2} > ~{sample}.~{j}.sam
    >>>

    output {
        File samFile = "~{sample}.~{j}.sam"
    }

    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        disks: "local-disk 100 SSD"
        memory: "8G"
        cpu: 2
    }
}

task Samblaster {
    input {
        String sample
        File samFile
    }

    command <<<
        /usr/local/bin/samblaster -a --addMateTags -i ~{samFile} -o ~{sample}.blast.sam
    >>>

    output {
        File blastsamFile = "~{sample}.blast.sam"
    }

    runtime {
        docker: "mgibio/alignment_helper-cwl:1.0.0"
        disks: "local-disk 100 SSD"
        memory: "8G"
        cpu: 2
    }

}

task toBam {
    input {
        String sample
        Int j
        File blastsamFile
    }

    command <<<
        /usr/local/bin/samtools sort -@ -8 -o ~{sample}.~{j}.final.bam ~{blastsamFile}
    >>>

    output {
        File bamFile = "~{sample}.~{j}.final.bam"
    }

    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        disks: "local-disk 100 SSD"
        memory: "8G"
        cpu: 2
    }
}
