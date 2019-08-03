# Mutect2: https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.4/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php
version 1.0

workflow runMutect2 {
    input {
        File interval_list
        File reference_fasta
        File gnomad_vcf
        File gnomad_vcf_index
        String tumor_sample
        String normal_sample
        File tumor_bam
        File normal_bam
    }

    Array[Array[String]] intervals = read_tsv(interval_list)

    scatter (interval in intervals) {
        Int i = interval[0]
        File interval_file = interval[1]

        call Mutect {
            input:
                interval_file = interval_file,
                reference_fasta = reference_fasta,
                gnomad_vcf = gnomad_vcf,
                gnomad_vcf_index = gnomad_vcf_index,
                tumor_sample = tumor_sample,
                normal_sample = normal_sample,
                tumor_bam = tumor_bam,
                normal_bam = normal_bam,
                i = i
        }
    
    }

    output {
        Array[File] vcfFile = Mutect.vcfFile
    }
}

task Mutect {
    input {
        File interval_file
        File reference_fasta
        File gnomad_vcf
        File gnomad_vcf_index
        String tumor_sample
        String normal_sample
        File tumor_bam
        File normal_bam
        Int i
    }

    command <<<
        /usr/gitc/gatk4/gatk-launch --javaOptions "-Xmx4g" Mutect2 \
          -R ~{reference_fasta} \
          -I ~{tumor_bam} \
          -tumor ~{tumor_sample} \
          -I ~{normal_bam} \
          -normal ~{normal_sample} \
          --germline_resource ~{gnomad_vcf} \
          -L ~{interval_file} \
          -O ~{tumor_sample}.~{i}.mutect.snv_indel.vcf.gz
    >>>

    output {
        File vcfFile = "~{tumor_sample}.~{i}.mutect.snv_indel.vcf.gz"
    }

    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        disks: "local-disk 100 SSD"
        memory: "16G"
        cpu: 1
    }
}