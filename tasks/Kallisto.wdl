version 1.0

workflow runKallisto {
    input {
        String sample
        Array[File] fastqList
        File kallisto_index
    }

    call runQuant {
        input:
            sample = sample,
            fastqList = fastqList,
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
        Array[File]+ fastqList
        File kallisto_index
    }

    command <<<
        /usr/local/bin/kallisto quant -i ~{kallisto_index} -b 100 --fusion --fr-stranded -o ~{sample}.kallisto ~{sep=" " fastqList}

        for i in ~{sample}/*; do new=`echo $i | tr '/' '.'`; mv $i $new; done
    >>>

    output {
        Array[File] kallistoOut = glob("~{sample}.kallisto*")
        File pizzlyInput = "~{sample}.fusion.txt"
    }

    runtime {
        docker: "zlskidmore/kallisto:latest"
        disks: "local-disk 100 SSD"
        memory: "8G"
        cpu: 2
    }
}
