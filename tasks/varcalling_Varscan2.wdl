# Varscan2: https://github.com/genome/analysis-workflows/blob/master/definitions/tools/varscan_somatic.cwl
version 1.0

task Varscan {
    input {
        File tumor_bam
        File tumor_bam_index
        String tumor_sample
        File normal_bam
        File normal_bam_index
        File reference_fasta
        File reference_fasta_index
    }

    command <<<
    java -jar /opt/varscan/VarScan.jar somatic \
    <(/opt/samtools/bin/samtools mpileup --no-baq -f ~{reference_fasta} ~{normal_bam} ~{tumor_bam}) \
    ~{tumor_sample} \
    --strand-filter 0 \
    --min-coverage 20 \
    --min-var-freq 0.05 \
    --p-value 0.99 \
    --mpileup 1 \
    --output-vcf 1
    >>>

    output {
        File snpFile = "~{tumor_sample}.snp.vcf"
        File indelFile = "~{tumor_sample}.indel.vcf"
    }

    runtime {
        docker: "mgibio/cle:v1.3.1"
        disks: "local-disk 100 SSD"
        memory: "16G"
        cpu: 1
    }
}