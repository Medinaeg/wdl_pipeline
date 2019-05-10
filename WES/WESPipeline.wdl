version 1.0

##TODO: Add Fastq2 NA if needed/does not exist
import "./tasks/runBWA.wdl" as BWA
import "./tasks/MergeAlignedBams.wdl" as MergeAlignedBams

#import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/runBWA.wdl" as BWA
#import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/MergeAlignedBams.wdl" as MergedAlignedBams

workflow myWorkflow {
    input {
        File fofn
        String reference_prefix
        File reference_fa
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

        call  BWA.runBWA as BWA {
            input:
                sample = sample,
                fileList = getSamplesPerIndex.pairedFileList,
                reference_prefix = reference_prefix,
                reference_fa = reference_fa
        }

        if ( getSamplesPerIndex.nPairsOfFastqs != "1" ) {
            call MergeAlignedBams.mergeBams as mergeBams {
                input:
                    sample = sample,
                    bamFiles = BWA.bamFile
            }
        }

       File ouputAlignedBam = select_first([mergeBams.mergedBam, BWA.bamFile])
    }

    output {
        Array[File] outputAlignedBams = ouputAlignedBam
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
        Int i
        String sample
        File fofn
    }

    command <<<
        awk -v s="~{sample}" 'BEGIN{OFS="\t";} $1 == s {print $0,NR}' ~{fofn} > STDOUT.~{i}
        cat STDOUT.~{i} | cut -f2-3 | tr '\t' '\n' > FILELIST.~{i}
        wc -l STDOUT.~{i} > NLINES.~{i}
    >>>

    output {
        Array[Array[String]] pairedFileList = read_tsv("STDOUT.~{i}")
        Array[File] fastqList = read_lines("FILELIST.~{i}")
        String nPairsOfFastqs = read_string("NLINES.~{i}")
    }
}