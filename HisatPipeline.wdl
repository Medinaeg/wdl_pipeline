version 1.0

##TODO: Add Fastq2 NA if needed/does not exist
#import "./tasks/HisatAlignment.wdl" as Hisat
#import "./tasks/MergeAlignedBams.wdl" as MergeAlignedBams

import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/HisatAlignment.wdl" as Hisat
import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/MergeAlignedBams.wdl" as MergeAlignedBams

workflow myWorkflow {
    input {
        File fofn
        String strandness

        String hisat_prefix
        File hisat_index
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

        call Hisat.runAlignments as Hisat {
            input:
                sample = sample,
                fileList = getSamplesPerIndex.pairedFileList,
                strandness = strandness,
                hisatPrefix = hisat_prefix,
                hisatIndex = hisat_index
        }

        if ( getSamplesPerIndex.nPairsOfFastqs > 1 ) {
            call MergeAlignedBams.mergeBams as mergeBams {
                input:
                    sample = sample,
                    bamFiles = Hisat.bamFile
            }
        }

        File outputAlignedBam = select_first([mergeBams.mergedBam, Hisat.bamFile[0]])
    }

    output {
        Array[File] outputAlignedBams = outputAlignedBam
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

task getSamplesPerIndex {
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
