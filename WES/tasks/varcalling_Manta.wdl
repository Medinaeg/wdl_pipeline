version 1.0

task runManta {
    input {
        File normal_bam
        File normal_bamindex
        File tumor_bam
        File tumor_bamindex
        File reference_fasta
        File reference_fasta_index
        String tumor_sample
    }

    command <<<
        dir=$(echo ~{tumor_bam} | sed 's/~{tumor_sample}.FINAL.bam//')
        mv ~{tumor_bamindex} $dir
        mv ~{normal_bamindex} $dir

        /manta-1.4.0.centos6_x86_64/bin/configManta.py --normalBam ~{normal_bam} --tumorBam ~{tumor_bam} --referenceFasta ~{reference_fasta} --exome --runDir $PWD
        $PWD/runWorkflow.py -m local -j 3 -g 10
        mv results/variants/diploidSV.vcf.gz results/variants/~{tumor_sample}.Manta.diploidSV.vcf.gz
        mv results/variants/somaticSV.vcf.gz results/variants/~{tumor_sample}.Manta.somaticSV.vcf.gz
        mv results/variants/candidateSV.vcf.gz results/variants/~{tumor_sample}.Manta.candidateSV.vcf.gz
        mv results/variants/candidateSmallIndels.vcf.gz results/variants/~{tumor_sample}.Manta.candidateSmallIndels.vcf.gz
    >>>

    output {
        File diploidFile = "results/variants/~{tumor_sample}.Manta.diploidSV.vcf.gz"
        File somaticFile = "results/variants/~{tumor_sample}.Manta.somaticSV.vcf.gz"
        File candidateFile = "results/variants/~{tumor_sample}.Manta.candidateSV.vcf.gz"
        File candidateSmallIndelsFile = "results/variants/~{tumor_sample}.Manta.candidateSmallIndels.vcf.gz"
    }

    runtime {
        docker: "kfdrc/manta"
        disks: "local-disk 600 SSD"
        memory: "32G"
        cpu: 4
    }
    
}
