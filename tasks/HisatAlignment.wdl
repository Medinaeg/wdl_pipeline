version 1.0

workflow runAlignments {
    input {
        String sample
        File fastq1
        File fastq2
        String strandness
        File hisatIndex
    }

    call getReadInfo {
        input:
            sample = sample,
            fastq1 = fastq1
    }

    call hisatCommand {
        input:
            sample = sample,
            fastq1 = fastq1,
            fastq2 = fastq2,
            hisatIndex = hisatIndex,
            strandness = strandness,
            id = getReadInfo.FastqInfo[0],
            pu = getReadInfo.FastqInfo[1],
            sm = getReadInfo.FastqInfo[2],

    }

    output {
        File bamFile = hisatCommand.bamFile
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
        String sample
        File fastq1
        File fastq2
        File hisatIndex
        String strandness
        String id
        String pu
        String sm
    }

    command <<<
        if [ "~{strandness}" == "NA" ]; then
            /usr/local/bin/hisat2 -x ~{hisatIndex} --rg-id ~{id} --rg PL:ILLUMINA --rg PU:~{sample} --rg LB:~{id}.~{sm} --rg SM:~{sample} -1 ~{fastq1} -2 ~{fastq2} -S ~{sample}.align.sam
        else
            /usr/local/bin/hisat2 -x ~{hisatIndex} --rg-id ~{id} --rg PL:ILLUMINA --rg PU:~{sample} --rg LB:~{id}.~{sm} --rg SM:~{sample} --rna-strandness ~{strandness} -1 ~{fastq1} -2 ~{fastq2} -S ~{sample}.align.sam
        fi

        samtools sort -@ -8 -o ~{sample}.final.bam ~{sample}.align.sam     
    >>>

     output{
         File bamFile =  "~{sample}.final.bam"
     }

    runtime {
        docker: "limesbonn/hisat2:latest"
        disks: "local-disk 100 SSD"
        memory: "8G"
        cpu: 2
    }
}
