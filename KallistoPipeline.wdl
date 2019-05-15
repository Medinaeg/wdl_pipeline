version 1.0

##TODO: Add Fastq2 NA if needed/does not exist
import "./tasks/Kallisto.wdl" as Kallisto
import "./tasks/Pizzly.wdl" as Pizzly

#import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/Kallisto.wdl" as Kallisto
#import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/Pizzly.wdl" as Pizzly

workflow myWorkflow {
    input {
        File fofn
        String strandness

        File reference_gtf
        File kallisto_index
        File reference_cdna
    }

    call splitSamples {
        input:
            fofn = fofn
        }

    scatter (sampleIndex in splitSamples.sampleList) {
        String index = sampleIndex[0]
        String sample = sampleIndex[1]

        call getSamplesPerIndex {
            input:
                i = index,
                sample = sample,
                fofn = fofn
        }

        call Kallisto.runKallisto as runKallisto {
            input:
                sample = sample,
                fastqList = getSamplesPerIndex.fastqList,
                kallisto_index = kallisto_index
        }

        call Pizzly.getFusions as runPizzly {
            input:
                sample = sample,
                reference_gtf = reference_gtf,
                reference_cdna = reference_cdna,
                pizzlyInput = runKallisto.pizzlyInput
        }
    }

    output {
        Array[Array[File]] kallistoOut = runKallisto.kallistoOut
        Array[Array[File]] pizzlyOut = runPizzly.pizzlyOut
    }
}

task splitSamples {
    input {
        File fofn
    }

    command <<<
        cut -f1 ~{fofn} | sort | uniq | awk 'BEGIN{OFS="\t";} {print NR,$0}' > STDOUT
    >>>

    output {
        Array[Array[String]] sampleList = read_tsv("STDOUT")
    }
}

task getFastqs {
    input {
        String i
        String sample
        File fofn
    }

    command <<<
        awk -v s="~{sample}" '$1 == s {print}' ~{fofn} > STDOUT.~{i}
        cat STDOUT.~{i} | cut -f2-3 | tr '\t' '\n' > FILELIST.~{i}
        wc -l STDOUT.~{i} > NLINES.~{i}
    >>>

    output {
        Array[Array[String]] pairedFileList = read_tsv("STDOUT.~{i}")
        Array[File] fastqList = read_lines("FILELIST.~{i}")
        Int nPairsOfFastqs = read_string("NLINES.~{i}")
    }
}
