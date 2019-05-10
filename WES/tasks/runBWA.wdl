version 1.0

workflow runBWA {
    input {
        String sample
        Array[Array[String]] fileList #tsv of sample, fastq1/2, index
        String reference_prefix
        File reference_fa
    }

    scatter (filePair in fileList) {
        String sample = filePair[0]
        File fastq1 = filePair[1]
        File fastq2 = filePair[2]
        Int j = filePair[3]

        call getReadInfo {
            input:
                sample = sample,
                fastq1 = fastq1
        }

        call BWACommand {
            input:
                j = j,
                sample = sample,
                fastq1 = fastq1,
                fastq2 = fastq2,
                reference_prefix = reference_prefix,
                reference_fa = reference_fa,
                id = getReadInfo.FastqInfo[0],
                pu = getReadInfo.FastqInfo[1],
                sm = getReadInfo.FastqInfo[2],
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

task getReadInfo {
    input {
        String sample
        File fastq1
    }

    command <<<
        zcat < ~{fastq1} | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/./g' >> FastqInfo
        zcat < ~{fastq1} | head -n 1 | cut -f 3-4 -d":" | sed 's/@//' | sed 's/:/./g' >> FastqInfo
        zcat < ~{fastq1} | head -n 1 | grep -Eo "[ATGCN]+$" >> FastqInfo
    >>>

    output {
        Array[String] FastqInfo = read_lines("FastqInfo")
    }

    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        disks: "local-disk 100 SSD"
        memory: "8G"
        cpu: 2
    }
}


task BWACommand {
    input {
        Int j
        String sample
        File fastq1
        File fastq2
        String reference_prefix
        File reference_fa
        String id
        String pu
        String sm  
    }

    String bwaDetails = 'echo "@RG\tID:~{id}\tPU:~{pu}"."~{sm}\tSM:~{sample}\tLB:~{id}"."~{sm}\tPL:ILLUMINA\tCN:UCLA"'
    
    command <<<
        bwa mem -K 100000000 -t 8 -R ~(bwaDetails) ./~{reference_prefix} ~{fastq1} ~{fastq2} > ~{sample}.~{j}.sam 
#        bwa mem -K 100000000 -t 8 -R $(echo "@RG\tID:$id\tPU:$pu"."$sm\tSM:$sample\tLB:$id"."$sm\tPL:ILLUMINA\tCN:UCLA") ${reference} ${r1file} ${r2file} > ${prefix}.sam
    >>>

    output {
        File samFile = "~{sample}.~{j}.sam"
    }

    runtime {
        docker: "zlskidmore/hisat2:latest"
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
        samblaster -a --addMateTags -i ~{samFile} -o ${sample}.blast.sam
        /u/flashscratch/k/katiecam/software/samblaster/samblaster -a --addMateTags -i ${prefix}.sam -o ${prefix}.blast.sam
    >>>

    output {
        File blastsamFile = "~{sample}.blast.sam"
    }


}

task toBam {
    input {
        String sample
        Int j
        File blastsamFile
    }

    command <<<
        samtools sort -@ -8 -o ~{sample}.~{j}.final.bam ~{blastsamFile}
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