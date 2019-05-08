version 1.0

##TODO: Add Fastq2 NA if needed/does not exist
import "./tasks/Kallisto.wdl" as Kallisto
import "./tasks/Pizzly.wdl" as Pizzly

#import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/Kallisto.wdl" as Kallisto
#import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/Pizzly.wdl" as Pizzly

workflow myWorkflow {
    input {
        String
        String strandness

        String hisat_index
        File reference_gtf
        File kallisto_index
        File reference_cdna
    }

    Array[Array[String]] inputSamples = read_tsv(fofn)

    scatter (line in inputSamples) {
        String sample = line[0]
        File fastq1 = line[1]
        File fastq2 = line[2]

        call Kallisto.runKallisto as runKallisto {
            input:
                sample = sample,
                fastq1 = fastq1,
                fastq2 = fastq2,
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
        Array[File] quantFile = runKallisto.quantFile
        Array[File] pizzlyInput = runKallisto.pizzlyInput
        Array[File] pizzlyOutput = runPizzly.unfilteredJSON
    }
}
