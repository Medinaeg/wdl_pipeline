version 1.0

workflow getFusions {
    input {
        String sample
        File reference_gtf
        File reference_cdna
        File pizzlyInput
    }

    call runPizzly {
        input:
            sample = sample,
            reference_gtf = reference_gtf,
            reference_cdna = reference_cdna,
            pizzlyInput = pizzlyInput
    }

    output {
        File unfilteredJSON = runPizzly.unfilteredJSON
    }
}

task runPizzly {
    input {
        String sample
        File reference_gtf
        File reference_cdna
        File pizzlyInput
    }

    command <<<
    /usr/local/bin/pizzly -k 31 --gtf ~{reference_gtf} --align-score 2 --insert-size 400 --fasta ~{reference_cdna} --output ~{sample}/pizzly ~{pizzlyInput}
    >>>

    output {
        File unfilteredJSON = "~{sample}/pizzly.unfiltered.json"
    }

    runtime {
        docker: "chrisamiller/docker-pizzly:latest"
        disks: "local-disk 100 SSD"
        memory: "8G"
        cpu: 2
    }
}
