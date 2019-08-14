version 1.0

# import "./tasks/varcalling_Varscan2_fullprocessing.wdl" as Varscan
# import "./tasks/varcalling_Strelka_fullprocessing.wdl" as Strelka
# import "./tasks/mergeVariantCallers.wdl" as MergeVCFs
# import "./tasks/varcalling_SomaticSniper.wdl" as SomaticSniper
# import "./tasks/createSequenza.wdl" as Sequenza
# import "./tasks/varcalling_Manta.wdl" as Manta

import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/WES/tasks/varcalling_Varscan2_revision.wdl" as Varscan
import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/WES/tasks/varcalling_Strelka_revision.wdl" as Strelka
import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/WES/tasks/mergeVariantCallers_snvs.wdl" as MergeSNVs
import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/WES/tasks/mergeVariantCallers_indels.wdl" as MergeINDELs
import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/WES/tasks/varcalling_Manta_revision.wdl" as Manta

workflow SomaticVaraintDetection {
    input {
        File fofn_bams_paired
        File pathsToReferenceFastaFiles
        File gcWiggle
    }
# Prelim step: Convert fofn_bams_paired to Array for scatter
    Array[Array[String]] map_bams = read_tsv(fofn_bams_paired)
# Prelim step: Get reference.fa files as bundle
    call getReferenceFiles {
        input:
            pathsToReferenceFastaFiles = pathsToReferenceFastaFiles
    }
#
    scatter (pair in map_bams) {
        String tumorSample = pair[0]
        File tumorBam = pair[1]
        String normalSample = pair[2]
        File normalBam = pair[3]
# 1. Gather information of bam through flagstat and indexes bams
        call processBam {
            input:
                tumorSample = tumorSample,
                tumorBam = tumorBam,
                normalSample = normalSample,
                normalBam = normalBam
        }

# 2. Run Varscan on paired tumor/normal
        call Varscan.runVarscan as Varscan {
            input:
                tumor_bam = tumorBam,
                tumor_bam_index = processBam.tumorBamIndex,
                tumor_sample = tumorSample,
                normal_bam = normalBam,
                normal_bam_index = processBam.normalBamIndex,
                referenceFastaFiles = getReferenceFiles.referenceFastaFiles
        }
# 3. Run Strelks on paired tumor/normal
        call Strelka.runStrelka as Strelka {
            input:
                normalbam = normalBam,
                normalbamindex = processBam.normalBamIndex,
                tumorsample = tumorSample,
                tumorbam = tumorBam,
                tumorbamindex = processBam.tumorBamIndex,
                referenceFastaFiles = getReferenceFiles.referenceFastaFiles
        }

        call MergeSNVs.runCombineVariantsSNVs as MergeSNVs {
            input:
                tumor_sample = tumorSample,
                varscanFile = Varscan.finalVarscanSnp,
                varscanFileIndex = Varscan.finalVarscanSnpindex,
                strelkaFile = Strelka.snvFile,
                strelkaFileIndex = Strelka.snvFileIndex

        }

        call MergeINDELs.runCombineVariantsINDELs as MergeINDELs {
            input:
                tumor_sample = tumorSample,
                varscanFile = Varscan.finalVarscanIndel,
                varscanFileIndex = Varscan.finalVarscanIndelindex,
                strelkaFile = Strelka.indelFile,
                strelkaFileIndex = Strelka.indelFileIndex          
        }

# 6. Run Manta on paired tumor/normal
        call Manta.runManta as Manta {
            input:
                normal_bam = normalBam,
                normal_bamindex = processBam.normalBamIndex,
                tumor_bam = tumorBam,
                tumor_bamindex = processBam.tumorBamIndex,
                tumor_sample = tumorSample,
                referenceFastaFiles = getReferenceFiles.referenceFastaFiles
        }

    }
# 7. Ouputs: flagstat information, bam index, Varscan vcf, Strelka vcf, Somatic Sniper vcf, Manta vcf
    output {
        Array[File] outputtumorFlagstat = processBam.tumorFlagstat
        Array[File] outputnormalFlagstat = processBam.normalFlagstat
        Array[File] outputtumorBamIndex = processBam.tumorBamIndex
        Array[File] outputnormalBamIndex = processBam.normalBamIndex 
        Array[Array[File]] outputVarscanFiles = Varscan.allVarscanFiles
        Array[File] outputStrelkasnvFiles = Strelka.snvFile
        Array[File] outputStrelkaindelFiles = Strelka.indelFile
        Array[File] outputSNVmergedFiles = MergeSNVs.mergevcf
        Array[File] outputSNVmergedFilesIndex = MergeSNVs.mergevcfindex
        Array[File] outputINDELmergedFiles = MergeINDELs.mergevcf
        Array[File] outputINDELmergedFilesIndex = MergeINDELs.mergevcfindex
        Array[Array[File]] outputMantaFiles = Manta.Mantaoutput
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

task processBam {
    input {
        String tumorSample 
        File tumorBam
        String normalSample
        File normalBam
    }

    command <<<
        /usr/local/bin/samtools flagstat ~{tumorBam} > ~{tumorSample}.flagstat.txt
        
        /usr/local/bin/samtools flagstat ~{normalBam} > ~{normalSample}.flagstat.txt

        /usr/local/bin/samtools index ~{tumorBam} $PWD/~{tumorSample}.FINAL.bam.bai

        /usr/local/bin/samtools index ~{normalBam} $PWD/~{normalSample}.FINAL.bam.bai
    >>>

    output {
        File tumorFlagstat = "~{tumorSample}.flagstat.txt"
        File normalFlagstat = "~{normalSample}.flagstat.txt"
        File tumorBamIndex = "~{tumorSample}.FINAL.bam.bai"
        File normalBamIndex = "~{normalSample}.FINAL.bam.bai"
    }

#    Remove comments if running locally
    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        disks: "local-disk 100 SSD"
        memory: "16G"
        cpu: 1
    }
}