version 1.0

##TODO: Add Fastq2 NA if needed/does not exist
import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/WES/tasks/PicardMarkDuplicatesBQSR.wdl" as MarkDuplicatesBQSR

workflow myWorkflow {
    input {
        File fofn_bams
        File pathsToReferenceFastaFiles
        File pathsToGATK4ReferenceFiles
    }
# Prelim step: Get reference.fa files as bundle
    call getReferenceFiles {
        input:
            pathsToReferenceFastaFiles = pathsToReferenceFastaFiles,
            pathsToGATK4ReferenceFiles = pathsToGATK4ReferenceFiles
    }

    Array[Array[String]] map_bams = read_tsv(fofn_bams)

    scatter (sampleBam in map_bams) {
        String sample = sampleBam[0]
        File bamFile = sampleBam[1]

       call MarkDuplicatesBQSR.GatkCommands as MDBQSR {
            input:
                sample = sample,
                bamFile = bamFile,
                referenceFastaFiles = getReferenceFiles.referenceFastaFiles,
                referenceGATK4Files = getReferenceFiles.referenceGATK4Files
       }

    }
# 2D. Output aligned BAM files
    output {
        Array[File] outputFinalBams = MDBQSR.finalBam
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
