version 1.0

workflow StringTieFPKM {
    input {
        String sample
        File alignedBam
        File reference_gtf
    }

    call StringTieCommand {
        input:
            sample = sample,
            alignedBam = alignedBam,
            gtf = reference_gtf
    }

    output {
        File fpkmFile = StringTieCommand.fpkmFile
    }
}


task StringTieCommand {
    input {
        String sample
        File alignedBam
        File gtf
    }

    command <<<
        /usr/local/bin/stringtie -G ~{gtf} -e -B -o ~{sample}.transcripts.gtf -A ~{sample}.abundances.tsv ~{alignedBam}
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
