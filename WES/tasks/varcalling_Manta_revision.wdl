version 1.0

task runManta {
    input {
        File normal_bam
        File normal_bamindex
        File tumor_bam
        File tumor_bamindex
        String tumor_sample
        Array[File] referenceFastaFiles
    }

    File reference_fasta = referenceFastaFiles[0]

    command <<<
        dir=$(echo ~{tumor_bam} | sed 's/~{tumor_sample}.FINAL.bam//')
        mv ~{tumor_bamindex} $dir
        mv ~{normal_bamindex} $dir

        /manta-1.4.0.centos6_x86_64/bin/configManta.py --normalBam ~{normal_bam} --tumorBam ~{tumor_bam} --referenceFasta ~{reference_fasta} --exome --runDir $PWD
        $PWD/runWorkflow.py -m local -j 3 -g 10
        mv results/variants/diploidSV.vcf.gz results/variants/~{tumor_sample}_Manta.diploidSV.vcf.gz
        mv results/variants/somaticSV.vcf.gz results/variants/~{tumor_sample}_Manta.somaticSV.vcf.gz
        mv results/variants/candidateSV.vcf.gz results/variants/~{tumor_sample}_Manta.candidateSV.vcf.gz
        mv results/variants/candidateSmallIndels.vcf.gz results/variants/~{tumor_sample}_Manta.candidateSmallIndels.vcf.gz
    >>>

    output {
        Array[File] Mantaoutput = glob("*vcf*")
    }

    runtime {
        docker: "kfdrc/manta"
        disks: "local-disk 600 SSD"
        memory: "32G"
        cpu: 4
    }
    
}
