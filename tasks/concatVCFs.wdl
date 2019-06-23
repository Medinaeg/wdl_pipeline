version 1.0

task concatVCFs {
    input{
        String tumor_sample
        Array[File] vcfFiles
    }

    command <<<
        vcf-concat ~{sep=" " vcfFiles} | gzip -c > ~{tumor_sample}.FINAL.mutect.snv_indel.vcf.gz
    >>>

    output {
        File concatVCFfiles = "~{tumor_sample}.FINAL.mutect.snv_indel.vcf.gz"
    }

    runtime {
        docker: "opengenomics/vcftools-tools:latest"
        disks: "local-disk 100 SSD"
        memory: "16G"
        cpu: 1  
    }
}