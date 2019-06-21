# Mutect2: https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.4/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php
version 1.0

workflow runAlignments {
    input {
        Array[File]+ interval_list
        File reference_fasta
        File tumor_bam
        String tumor_sample
        File normal_bam
        String normal_sample
        File gnomad_vcf
        File interval
        Int i
    }

    scatter (interval in interval_list) {
        Int i = interval[0]
        File interval_file = interval[1]

        call Mutect {
            input:
                reference_fasta = reference_fasta,
                tumor_bam = tumor_bam,
                tumor_sample = tumor_sample,
                normal_bam = normal_bam,
                normal_sample = normal_sample,
                gnomad_vcf = gnomad_vcf,
                interval = interval_file,
                i = i
        }
    }

    output {

    }
}

task Mutect {
    input {
        File reference_fasta
        File tumor_bam
        String tumor_sample
        File normal_bam
        String normal_sample
        File gnomad_vcf
        File interval
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
          -L ~{interval} \
          -O ~{tumor_sample}.~{i}.mutect.snv_indel.vcf.gz
    >>>

    output {
        File "~{tumor_sample}.~{i}.mutect.snv_indel.vcf.gz"
    }

    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        disks: "local-disk 100 SSD"
        memory: "16G"
        cpu: 1
    }
}

task aggMutect {
    input {
        Array[File]+
    }
}
