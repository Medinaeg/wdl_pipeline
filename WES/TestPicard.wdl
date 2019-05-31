version 1.0

##TODO: Add Fastq2 NA if needed/does not exist
#import "./tasks/PicardMarkDuplicatesBQSR.wdl" as MarkDuplicatesBQSR

import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/WES/tasks/PicardMarkDuplicatesBQSR.wdl" as MarkDuplicatesBQSR

workflow myWorkflow {
    input {
        File bamFile
        String sample
        File pathsToReferenceFastaFiles
        File pathsToGATK4ReferenceFiles
    }

    call getReferenceFiles {
        input:
            pathsToReferenceFastaFiles = pathsToReferenceFastaFiles,
            pathsToGATK4ReferenceFiles = pathsToGATK4ReferenceFiles
    }

    call MarkDuplicatesBQSR.GatkCommands as MDBQSR {
            input:
                sample = sample,
                bamFile = bamFile,
                referenceFastaFiles = getReferenceFiles.referenceFastaFiles,
                referenceGATK4Files = getReferenceFiles.referenceGATK4Files
    }

    output {
        File outputFinalBams = MDBQSR.finalBam
    }
}

task getReferenceFiles {
    input {
        File pathsToReferenceFastaFiles
        File pathsToGATK4ReferenceFiles
    }

    command <<<
        cut -f1 ~{pathsToReferenceFastaFiles} >> STDOUT
        cut -f1 ~{pathsToGATK4ReferenceFiles} >> STDOUT1
    >>>

    output {
        Array[File] referenceFastaFiles = read_lines("STDOUT")
        Array[File] referenceGATK4Files = read_lines("STDOUT1")
    }
}