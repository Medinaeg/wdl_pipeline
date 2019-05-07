version 1.0

##TODO: Add Fastq2 NA if needed/does not exist
#import "./tasks/HisatAlignment.wdl" as Hisat
#import "./tasks/SamToBam.wdl" as SamToBam
#import "./tasks/HTSeq2version2.wdl" as HTseq
#import "./tasks/Stringtie.wdl" as Stringtie
#import "./tasks/Kallisto.wdl" as Kallisto
#import "./tasks/Pizzly.wdl" as Pizzly

import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/HisatAlignment.wdl" as Hisat
import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/SamToBam.wdl" as SamToBam
import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/HTSeq2version2.wdl" as HTseq
import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/Stringtie.wdl" as Stringtie
import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/Kallisto.wdl" as Kallisto
import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/Pizzly.wdl" as Pizzly

workflow myWorkflow {
    input {
        File fofn
        String outdir
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

        call Hisat.runAlignments as Hisat {
            input:
                outdir = outdir,
                sample = sample,
                fastq1 = fastq1,
                fastq2 = fastq2,
                strandness = strandness,
                hisatIndex = hisat_index
        }

        call SamToBam.SamToBAM as ToBam {
            input:
                sample = sample,
                outdir = outdir,
                samFile = Hisat.samFile
        }

        call HTseq.HTSeq2 as ToCounts {
            input:
                sample = sample,
                outdir = outdir,
                alignedBam = ToBam.finalBam,
                reference_gtf = reference_gtf
        }

        call Stringtie.StringTie as getFPKM {
            input:
                sample = sample,
                outdir = outdir,
                alignedBam = ToBam.finalBam,
                reference_gtf = reference_gtf
        }

        call Kallisto.runKallisto as runKallisto {
            input:
                sample = sample,
                outdir = outdir,
                fastq1 = fastq1,
                fastq2 = fastq2,
                kallisto_index = kallisto_index
        }

        call Pizzly.getFusions as runPizzly {
            input:
                sample = sample,
                outdir = outdir,
                reference_gtf = reference_gtf,
                reference_cdna = reference_cdna,
                pizzlyInput = runKallisto.pizzlyInput
        }
    }

    output {
        Array[File] alignedBam = ToBam.finalBam
        Array[File] countsFile = ToCounts.countsFile
        Array[File] fpkmFile = getFPKM.fpkmFile
        Array[File] quantFile = runKallisto.quantFile
        Array[File] pizzlyInput = runKallisto.pizzlyInput
        Array[File] pizzlyOutput = runPizzly.unfilteredJSON
    }
}
