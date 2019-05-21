version 1.0

workflow HTSeq2 {
    input {
        String sample
        File alignedBam
        File reference_gtf
    }

    call Counts {
        input:
            sample = sample,
            alignedBam = alignedBam,
            gtf = reference_gtf
    }

    output {
        File countsFile = Counts.countsFile
    }
}

task Counts {
    input {
        String sample
        File alignedBam
        File gtf
    }

    command <<<
        /usr/local/HTSeq-0.6.1p1/scripts/htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id ~{alignedBam} ~{gtf} > ~{sample}.HTSeq2counts.tsv
    >>>

    output {
        File countsFile = "~{sample}.HTSeq2counts.tsv"
    }

    runtime {
        docker: "dmccloskey/htseq-count:latest"
        disks: "local-disk 100 SSD"
        memory: "8G"
        cpu: 2
    }
}
