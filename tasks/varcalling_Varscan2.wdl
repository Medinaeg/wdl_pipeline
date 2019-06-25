# Varscan2: https://github.com/genome/analysis-workflows/blob/master/definitions/tools/varscan_somatic.cwl
version 1.0

###need to index BAM files
workflow runVarscan2 {
    input {
        File reference_fasta
        File tumor_bam
        File tumor_bam_index
        String tumor_sample
        File normal_bam
        File normal_bam_index
    }

    call Varscan {
        input:
            reference_fasta = reference_fasta,
            tumor_bam = tumor_bam,
            tumor_bam_index = tumor_bam_index,
            tumor_sample = tumor_sample,
            normal_bam = normal_bam,
            normal_bam_index = normal_bam_index
    }

    output {
        File snpFile = "~{tumor_sample}.snp.vcf"
        File indelFile = "~{tumor_sample}.indel.vcf"
    }
}

task Varscan {
    input {
        File tumor_bam
        File tumor_bam_index
        String tumor_sample
        File normal_bam
        File normal_bam_index
        File reference_fasta
    }

    #Change for specific need. These are default values.
    Int strand_filter = 0
    Int min_coverage = 20
    Float min_var_frequency = 0.05
    Float p_value = 0.99

    command <<<
        /usr/bin/varscan_helper.sh ~{tumor_bam} ~{normal_bam} ~{reference_fasta} ~{strand_filter} ~{min_coverage} ~{min_var_frequency} ~{p_value}
    >>>

    output {
        File snpFile = "~{tumor_sample}.snp.vcf"
        File indelFile = "~{tumor_sample}.indel.vcf"
    }

    runtime {
        docker: "mgibio/cle:v1.3.1"
        disks: "local-disk 300 SSD"
        memory: "16G"
        cpu: 1
    }
}