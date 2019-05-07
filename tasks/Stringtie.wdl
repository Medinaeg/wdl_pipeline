version 1.0

workflow StringTie {
    input {
        String sample
        String outdir
        File alignedBam
        File reference_gtf
    }

    call StringTie {
        input:
            sample = sample,
            outdir = outdir,
            bamFile = alignedBam,
            gtf = reference_gtf
    }

    output {
        File fpkmFile = StringTie.fpkmFile
    }
}


task StringTie{
    input {
        String sample
        String outdir
        File bamFile
        File gtf
    }

    command <<<
        /usr/local/bin/stringtie -G ~{gtf} -e -B -o ~{outdir}/expression/~{sample}.transcripts.gtf -A ~{outdir}/expression/~{sample}.abundances.tsv ~{bamFile}
    >>>

    output {
        File fpkmFile = "~{outdir}/expression/~{sample}.abundances.tsv"
    }

    runtime {
        docker: "ahwagner/stringtie:latest"
    }
}
