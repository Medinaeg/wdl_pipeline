version 1.0

workflow runKallisto {
    input {
        String sample
        String outdir
        File fastq1
        File fastq2
        File kallisto_index
    }

    call runQuant {
        input:
            sample = sample,
            outdir = outdir,
            fastq1 = fastq1,
            fastq2 = fastq2,
            kallisto_index = kallisto_index
    }

    output {
        File quantFile = runQuant.quantFile
        File pizzlyInput = runQuant.pizzlyInput
    }

}

task runQuant {
    input {
        String sample
        String outdir
        File fastq1
        File fastq2
        File kallisto_index
    }

    command <<<
    /usr/local/bin/kallisto quant -i ~{kallisto_index} -b 100 --fusion -o ~{sample} ~{fastq1} ~{fastq2}
    >>>

    output {
        File quantFile = "~{sample}/abundances.tsv"
        File pizzlyInput = "~{sample}/fusion.txt"
    }

    runtime {
        docker: "zlskidmore/kallisto:latest"
        disks: "local-disk 100 SSD"
        memory: "8G"
        cpu: 2
    }
}
