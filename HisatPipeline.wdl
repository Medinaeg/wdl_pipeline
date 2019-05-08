version 1.0

##TODO: Add Fastq2 NA if needed/does not exist
#import "./tasks/HisatAlignment.wdl" as Hisat
#import "./tasks/SamToBam.wdl" as SamToBam
#import "./tasks/HTSeq2version2.wdl" as HTseq
#import "./tasks/Stringtie.wdl" as Stringtie
#import "./tasks/Kallisto.wdl" as Kallisto
#import "./tasks/Pizzly.wdl" as Pizzly

import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/HisatAlignment.wdl" as Hisat

workflow myWorkflow {
    input {
        File fofn
        String outdir
        String strandness

        String hisat_prefix
        File hisat_index
        File reference_gtf
        File kallisto_index
        File reference_cdna
    }

    Array[Array[String]] inputSamples = read_tsv(fofn)

    scatter (line in inputSamples) {
        String sample = line[0]
        File fastq1 = line[1]
        File fastq2 = line[2]

        call Hisat.runAlignments as Hisat {
            input:
                sample = sample,
                fastq1 = fastq1,
                fastq2 = fastq2,
                strandness = strandness,
                hisatPrefix = hisat_prefix,
                hisatIndex = hisat_index
        }
    }

    output {
        Array[File] alignedBam = Hisat.bamFile
    }
}
