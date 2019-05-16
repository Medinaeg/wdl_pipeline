version 1.0

##TODO: Add Fastq2 NA if needed/does not exist

#import "./tasks/Kallisto.wdl" as Kallisto
#import "./tasks/Pizzly.wdl" as Pizzly
#import "./tasks/HisatAlignment.wdl" as Hisat
#import "./tasks/MergeAlignedBams.wdl" as MergeAlignedBams

import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/Kallisto.wdl" as Kallisto
import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/Pizzly.wdl" as Pizzly
import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/HisatAlignment.wdl" as Hisat
import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/MergeAlignedBams.wdl" as MergeAlignedBams


workflow myWorkflow {
    input {
        File fofn
        String strandness

        File reference_gtf
        File kallisto_index
        File reference_cdna
    }

# 1. Read in FOFN to get list of samples
    call listSamples {
        input fofn = fofn
    }
# Output: sampleList (TSV of index; sample)

# 2. Run pipeline on each Sample
    scatter( sampleIndex in listSamples.sampleList ) {
        String index = sampleIndex[0]
        String sample = sampleIndex[1]

        # A. Get fastq files per sample
        call listFastqFiles {
            input:
                i = index, # index of Sample
                sample = sample,
                fofn = fofn
        }
        # Output: pairedFileList (TSV of sample, Fastq1, Fastq2, pair Index)
        #       fastqList (Array of all Fastq files in order of [FQ1 FQ2 FQ1 FQ2])
        #       nPairsOfFastqs (Int; number of fastq pairs for the sample)

        # B. Run Kallisto on sample
        call Kallisto.runKallisto as runKallisto {
            input:
                sample = sample,
                fastqList = listFastqFiles.fastqList,
                kallisto_index = kallisto_index
        }
        # Output: kallistoOut (Array[File]), pizzlyInput (File)

        # C. Run Pizzly on Kallisto output
        call Pizzly.getFusions as runPizzly {
            input:
                sample = sample,
                reference_gtf = reference_gtf,
                reference_cdna = reference_cdna,
                pizzlyInput = runKallisto.pizzlyInput
        }
        # Output: pizzlyOut (Array[File])

        D. Run Hisat2 on sample (takes in TSV of Fastq Pairs)
        call Hisat.runAlignments as Hisat {
            input:
                sample = sample,
                fileList = listFastqFiles.pairedFileList,
                strandness = strandness,
                hisatPrefix = hisat_prefix,
                hisatIndex = hisat_index
        }
        # Output: bamFile (Array[File])
        # If there are >1 bam for the sample, merge:
        if ( listFastqFiles.nPairsOfFastqs > 1 ) {
            call MergeAlignedBams as mergeBams {
                input:
                    sample = sample,
                    bamFiles = Hisat.bamFile
            }
        }
        # If a merged bam was created, select. If not, pull the first element of Hisat.bamFile (since there is only a single element)
        File outputAlignedBam = select_first([mergeBams.mergedBam, Hisat.bamFile[0]])

        call renameFinalBam {
            input:
                sample = sample,
                bamFile = outputAlignedBam
        }
    }

    output {
        Array[Array[File]] kallistoResult = runKallisto.kallistoOut
        Array[Array[File]] pizzlyResult = runPizzly.pizzlyOut
        Array[File] finalBam = renameFinalBam.finalBam
    }
}

task listSamples {
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

task listFastqFiles {
    input {
        String i
        String sample
        File fofn
    }

    command <<<
        awk -v s="~{sample}" 'BEGIN{OFS="\t";} $1 == s {print $0,NR}' ~{fofn} > STDOUT.~{i}
        cat STDOUT.~{i} | cut -f2-3 | tr '\t' '\n' > FILELIST.~{i}
        cat STDOUT.~{i} | wc -l | sed 's/ //g' > NLINES.~{i}
    >>>

    output {
        Array[Array[String]] pairedFileList = read_tsv("STDOUT.~{i}")
        Array[File] fastqList = read_lines("FILELIST.~{i}")
        Int nPairsOfFastqs = read_int("NLINES.~{i}")
    }
}

task renameFinalBam {
    input {
        String sample
        File bamFile
    }

    command <<<
        mv bamFile ~{sample}.FINAL_ALIGNED.bam
    >>>

    output {
        File finalBam = "~{sample}.FINAL_ALIGNED.bam"
    }
}
