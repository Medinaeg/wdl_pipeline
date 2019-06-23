# Varscan2: https://github.com/genome/analysis-workflows/blob/master/definitions/tools/varscan_somatic.cwl
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
        File normal_bam
        Int strand_filter
        Float min_coverage
        Float min_var_frequency
        Float p_value
        File roi_bed
        Int i
    }

    command <<<
        ./usr/bin/varscan_helper.sh \
         ~{tumor_bam} \
         ~{normal_bam} \
         ~{reference_fasta} \
         ~{strand_filter} \
         ~{min_coverage} \
         ~{min_var_frequency} \
         ~{p_value} \
         ~{roi_bed}
    >>>

    output {
        File "~{tumor_sample}.~{i}.mutect.snv_indel.vcf.gz"
    }

    runtime {
        docker: ""
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
