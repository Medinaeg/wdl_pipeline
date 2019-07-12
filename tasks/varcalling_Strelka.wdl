version 1.0

task runStrelka {
    input {
        File normalBam
        File normalBamIndex
        String tumorSample
        File tumorBam
        File tumorBamIndex
        File referenceFasta
        File referenceFastafai
    }

#    File referenceFasta = referenceFastaFiles[0]

    command <<<
        /opt/strelka/bin/configureStrelkaSomaticWorkflow.py --normalBam=~{normalBam} --tumorBam=~{tumorBam} --referenceFasta=~{referenceFasta} --exome --runDir=$PWD
        python2 $PWD/runWorkflow.py -m local -j 3 -g 10
        mv results/variants/somatic.indels.vcf.gz results/variants/~{tumorSample}.Strelka.somatic.indels.vcf.gz
        mv results/variants/somatic.snvs.vcf.gz results/variants/~{tumorSample}.Strelka.somatic.snvs.vcf.gz        
    >>>

    output {
        File indelFile = "results/variants/~{tumorSample}.Strelka.somatic.indels.vcf.gz"
        File snvFile = "results/variants/~{tumorSample}.Strelka.somatic.snvs.vcf.gz"
    }

    runtime {
        docker: "mgibio/strelka-cwl:2.9.9"
        disks: "local-disk 400 SSD"
        memory: "32G"
        cpu: 4
    }
}