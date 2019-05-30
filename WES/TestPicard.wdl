version 1.0

##TODO: Add Fastq2 NA if needed/does not exist
#import "./tasks/PicardMarkDuplicatesBQSR.wdl" as MarkDuplicatesBQSR

import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/WES/tasks/PicardMarkDuplicatesBQSR.wdl" as MarkDuplicatesBQSR

workflow myWorkflow {
    input {
        File bamFile
        String sample
        File pathsToReferenceFastaFiles
        File thousG
        File thousGIndex
        File knownIndels
        File knownIndelsIndex
        File dbsnp
        File dbsnpIndex
    }

    call getReferenceFiles {
        input:
            pathsToReferenceFastaFiles = pathsToReferenceFastaFiles
    }

    call MarkDuplicatesBQSR.GatkCommands as MDBQSR {
            input:
                sample = sample,
                bamFile = bamFile,
                referenceFastaFiles = getReferenceFiles.referenceFastaFiles,
                thousG = thousG,
                knownIndels = knownIndels,
                dbsnp = dbsnp,
                thousGIndex = thousGIndex,
                knownIndelsIndex = knownIndelsIndex,
                dbsnp = dbsnp
    }

    output {
        File outputFinalBams = MDBQSR.finalBam
    }
}

task getReferenceFiles {
    input {
        File pathsToReferenceFastaFiles
    }

    command <<<
        cut -f1 ~{pathsToReferenceFastaFiles} >> STDOUT
    >>>

    output {
        Array[File] referenceFastaFiles = read_lines("STDOUT")
    }
}