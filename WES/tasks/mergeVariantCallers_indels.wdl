version 1.0

task runCombineVariantsINDELs {
    input {
        String tumor_sample
        File varscanFile
        File varscanFileIndex
        File strelkaFile
        File strelkaFileIndex
        Array[File] referenceFastaFiles
    }
    
    File referenceFasta = referenceFastaFiles[0]

    command <<<
        /home/biodocker/bin/GenomeAnalysisTK -T CombineVariants -R ~{referenceFasta} -genotypeMergeOptions PRIORITIZE --rod_priority_list varscan,strelka --variant:varscan ~{varscanFile} --variant:strelka ~{strelkaFile} -o ~{tumor_sample}_Varscan_Strelka_Merged_indels.FINAL.vcf
    >>>

    output {
        File mergevcf = "~{tumor_sample}_Varscan_Strelka_Merged_indels.FINAL.vcf"
        File mergevcfindex = "~{tumor_sample}_Varscan_Strelka_Merged_indels.FINAL.vcf.tbi"
    }

    runtime {
        docker: "nderoo324/gatk"
        disks: "local-disk 200 SSD"
        memory: "16G"
        cpu: 1
    }
}