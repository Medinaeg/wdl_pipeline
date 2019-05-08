version 1.0

workflow SamToBAM{
    input{
        String sample
        File samFile
        String outdir
    }
    call ToBam {
        input:
            sample = sample,
            samFile = samFile,
            outdir = outdir
    }
    call SortBam {
        input:
            sample = sample,
            alignedBam = ToBam.alignedBam,
            outdir = outdir
    }
    output {
        File finalBam = SortBam.finalBam
    }
}

task ToBam {
    input {
        String sample
        File samFile
        String outdir
    }

    command <<<
        samtools view -Sb -o ~{sample}.align.bam ~{samFile}
    >>>

    output {
        File alignedBam = "~{sample}.align.bam"
    }

    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
    }
}

task SortBam {
    input {
        String sample
        File alignedBam
        String outdir
    }

    command <<<
        samtools sort -o ~{sample}.final.bam ~{alignedBam}
    >>>

    output {
        File finalBam = "~{sample}.final.bam"
    }

    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        disks: "local-disk 100 SSD"
        memory: "8G"
        cpu: 2
    }
}
