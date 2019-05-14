version 1.0

##TODO: Add Fastq2 NA if needed/does not exist
import "./tasks/runBWA.wdl" as BWA
import "./tasks/MergeAlignedBams.wdl" as MergeAlignedBams
import "./tasks/PicardMarkDuplicatesBQSR.wdl" as MarkDuplicatesBQSR

#import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/runBWA.wdl" as BWA
#import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/MergeAlignedBams.wdl" as MergedAlignedBams
#import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/PicardMarkDuplicatesBQSR.wdl" as MarkDuplicatesBQSR

workflow myWorkflow {
    input {
        File fofn
        File reference_fa
        File thousG
        File knownIndels
        File dbsnp
    }
# 1. Get unique samples from file list (column 1 in fofn)
    call splitSamples {
        input:
            fofn = fofn
        }
# 2. Scatter across samples from sample list
    scatter (sampleIndex in splitSamples.sampleList) {
        String index = sampleIndex[0]
        String sample = sampleIndex[1]
# 2A. Get all pairs of fastq files per sample
        call getSamplesPerIndex {
            input:
                i = index,
                sample = sample,
                fofn = fofn
        }
# 2B. Run BWA-MEM alignment on pairs of FASTQs
        call  BWA.runBWA as BWA {
            input:
                sample = sample,
                fileList = getSamplesPerIndex.pairedFileList,
                reference_fa = reference_fa
        }
# 2C. Merge if there are more than one pair of FASTQ files
        if ( getSamplesPerIndex.nPairsOfFastqs != "1" ) {
            call MergeAlignedBams.mergeBams as mergeBams {
                input:
                    sample = sample,
                    bamFiles = BWA.bamFile
            }
        }

       File outputAlignedBam = select_first([mergeBams.mergedBam, BWA.bamFile])

       call MarkDuplicatesBQSR.GatkCommands as MDBQSR {
            input:
                sample = sample,
                bamFile = outputAlignedBam,
                reference_fa = reference_fa,
                thousG = thousG,
                knownIndels = knownIndels,
                dbsnp = dbsnp
       }

    }
# 2D. Output aligned BAM files
    output {
        Array[File] outputFinalBams = MDBQSR.finalBam
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
