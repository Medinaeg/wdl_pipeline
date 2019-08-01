version 1.0

task runSomaticSniper {
    input {
        File referenceFasta
        File referenceFastaIndex
        String tumor_Sample
        File tumor_Bam
        File normal_Bam
    }

    command <<<
        /usr/local/bin/bam-somaticsniper \
        -q 20 \
        -Q 20 \
        -F vcf \
        -f ~{referenceFasta} \
        ~{tumor_Bam} \
        ~{normal_Bam} \
        ~{tumor_Sample}_somaticsniper_variants.vcf
    >>>

    output {
        File somaticsniperVariants = "~{tumor_Sample}_somaticsniper_variants.vcf"
    }

    runtime {
        docker: "opengenomics/somatic-sniper"
        disks: "local-disk 200 SSD"
        memory: "16G"
        cpu: 1
    }       

}