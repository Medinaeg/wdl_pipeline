version 1.0

workflow runAlignments {
    input {
        String outdir
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
            outdir = outdir,
            sample = sample,
            files = if fastq2=="NA" then "-U fastq1" else "-1 fastq1 -2 fastq2",
            hisatIndex = hisatIndex,
            strandness = strandness,
            id = getReadInfo.FastqInfo[0],
            pu = getReadInfo.FastqInfo[1],
            sm = getReadInfo.FastqInfo[2],

    }

    output {
        File samFile = hisatCommand.samFile
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
        String outdir
        String sample
        String files
        File hisatIndex
        String strandness
        String id
        String pu
        String sm
    }

    command <<<
        if [ ~{strandness} == "NA" ]; then
            /usr/local/bin/hisat2 -x ~{hisatIndex} --rg-id ~{id} --rg PL:ILLUMINA --rg PU:~{sample} --rg LB:~{id}.~{sm} --rg SM:~{sample} ~{files} -S ~{outdir}/alignment/~{sample}.align.sam
        else
            /usr/local/bin/hisat2 -x ~{hisatIndex} --rg-id ~{id} --rg PL:ILLUMINA --rg PU:~{sample} --rg LB:~{id}.~{sm} --rg SM:~{sample} --rna-strandness ~{strandness} ~{files} -S ~{outdir}/alignment/~{sample}.align.sam
        fi
    >>>

     output{
         File samFile =  "~{outdir}/alignment/~{sample}.align.sam"
     }

    runtime {
        docker: "limesbonn/hisat2:latest"
    }
}
