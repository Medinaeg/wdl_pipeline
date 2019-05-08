version 1.0

workflow runAlignments {
    input {
        String sample
        Array[Array[String]] fileList # tsv of sample, fastq1/2, index
        String strandness
        String hisatPrefix
        File hisatIndex
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

        call hisatCommand {
            input:
                j = j,
                sample = sample,
                fastq1 = fastq1,
                fastq2 = fastq2,
                hisatPrefix = hisatPrefix,
                hisatIndex = hisatIndex,
                strandness = strandness,
                id = getReadInfo.FastqInfo[0],
                pu = getReadInfo.FastqInfo[1],
                sm = getReadInfo.FastqInfo[2],

        }
    }

    output {
        Array[File] bamFile = hisatCommand.bamFile
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
}

task hisatCommand {
    input {
        Int j
        String sample
        File fastq1
        File fastq2
        String hisatPrefix
        File hisatIndex
        String strandness
        String id
        String pu
        String sm
    }

    command <<<
        # If reference-data1/myhisat2.tar.gz is used, the directory will look like ./hisat2/GRCh38_HISAT
        tar -zxvf ~{hisatIndex} -C .
        # Set hisat_prefix to 'hisat2/GRCh38_HISAT2'

        if [ "~{strandness}" == "NA" ]; then
            /usr/local/bin/hisat2 -x ./~{hisatPrefix} --rg-id ~{id} --rg PL:ILLUMINA --rg PU:~{sample} --rg LB:~{id}.~{sm} --rg SM:~{sample} -1 ~{fastq1} -2 ~{fastq2} -S ~{sample}.align.sam
        else
            /usr/local/bin/hisat2 -x ./~{hisatPrefix} --rg-id ~{id} --rg PL:ILLUMINA --rg PU:~{sample} --rg LB:~{id}.~{sm} --rg SM:~{sample} --rna-strandness ~{strandness} -1 ~{fastq1} -2 ~{fastq2} -S ~{sample}.align.sam
        fi

        samtools sort -@ -8 -o ~{sample}.~{j}.final.bam ~{sample}.align.sam
    >>>

     output{
         File bamFile =  "~{sample}.~{j}.final.bam"
     }

    runtime {
        docker: "limesbonn/hisat2:latest"
        disks: "local-disk 100 SSD"
        memory: "8G"
        cpu: 2
    }
}
