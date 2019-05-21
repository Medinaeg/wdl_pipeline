version 1.0

workflow StringTie {
    input {
        String sample
        File alignedBam
        File reference_gtf
    }

    call StringTie {
        input:
            sample = sample,
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
        File bamFile
        File gtf
    }

    command <<<
        /usr/local/bin/stringtie -G ~{gtf} -e -B -o ~{sample}.transcripts.gtf -A ~{sample}.abundances.tsv ~{bamFile}
    >>>

    output {
        File fpkmFile = "~{sample}.abundances.tsv"
    }

    runtime {
        docker: "ahwagner/stringtie:latest"
        disks: "local-disk 100 SSD"
        memory: "8G"
        cpu: 2
    }
}
