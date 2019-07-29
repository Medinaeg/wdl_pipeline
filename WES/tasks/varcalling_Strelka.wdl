version 1.0

task runStrelka {
    input {
        File normalbam
        File normalbamindex
        String tumorsample
        File tumorbam
        File tumorbamindex
        File referenceFasta
        File referenceFastafai
    }

    command <<<
        dir=$(echo ~{tumorbam} | sed 's/~{tumorsample}.FINAL.bam//')
        mv ~{tumorbamindex} $dir
        mv ~{normalbamindex} $dir

        /opt/strelka/bin/configureStrelkaSomaticWorkflow.py --normalBam=~{normalbam} --tumorBam=~{tumorbam} --referenceFasta=~{referenceFasta} --exome --runDir=$PWD
        python2 $PWD/runWorkflow.py -m local -j 3 -g 10
        mv results/variants/somatic.indels.vcf.gz results/variants/~{tumorsample}.Strelka.somatic.indels.vcf.gz
        mv results/variants/somatic.snvs.vcf.gz results/variants/~{tumorsample}.Strelka.somatic.snvs.vcf.gz        
    >>>

    output {
        File indelFile = "results/variants/~{tumorsample}.Strelka.somatic.indels.vcf.gz"
        File snvFile = "results/variants/~{tumorsample}.Strelka.somatic.snvs.vcf.gz"
    }

    runtime {
        docker: "mgibio/strelka-cwl:2.9.9"
        disks: "local-disk 400 SSD"
        memory: "32G"
        cpu: 4
    }
}