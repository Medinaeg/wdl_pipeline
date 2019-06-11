version 1.0

workflow runAlignments {
    input {
        String sample
        Array[Array[String]] fileList # tsv of sample, fastq1/2, index
        String strandness
        String hisatPrefix
        Array[File] hisatIndex
    }

    scatter (filePair in fileList) {
        File fastq1 = filePair[1]
        File fastq2 = filePair[2]
        Int j = filePair[3]

        call runHisat {
            input:
                j = j,
                sample = sample,
                fastq1 = fastq1,
                fastq2 = fastq2,
                hisatPrefix = hisatPrefix,
                hisatIndex = hisatIndex,
                strandness = strandness
        }

        call toBam {
            input:
                j = j,
                sample = sample,
                samFile = runHisat.samFile
        }
    }

    output {
        Array[File] bamFile = toBam.bamFile
    }
}

task runHisat {
    input {
        Int j
        String sample
        File fastq1
        File fastq2
        String hisatPrefix
        Array[File]+ hisatIndex
        String strandness
    }

    command <<<
        id=$(zcat < ~{fastq1} | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/./g')
        pu=$(zcat < ~{fastq1} | head -n 1 | cut -f 3-4 -d":" | sed 's/@//' | sed 's/:/./g')
        sm=$(zcat < ~{fastq1} | head -n 1 | grep -Eo "[ATGCN]+$")

        /usr/local/bin/hisat2 -p 8 --dta -x ~{hisatPrefix} --rg-id $id --rg PL:ILLUMINA --rg PU:~{sample} --rg LB:$id --rg SM:~{sample} --rna-strandness ~{strandness} -1 ~{fastq1} ~{fastq2} -S ~{sample}.~{j}.align.sam
    >>>

    output {
        File samFile = "~{sample}.~{j}.align.sam"
    }

    runtime {
        docker: "zlskidmore/hisat2:latest"
        disks: "local-disk 100 SSD"
        memory: "16G"
        cpu: 1
    }
}

task toBam {
    input {
        String sample
        Int j
        File samFile
    }

    command <<<
        /usr/local/bin/samtools sort -@ -8 -o ~{sample}.~{j}.final.bam ~{samFile}
    >>>

    output {
        File bamFile = "~{sample}.~{j}.final.bam"
    }

    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        disks: "local-disk 100 SSD"
        memory: "8G"
        cpu: 1
    }
}
