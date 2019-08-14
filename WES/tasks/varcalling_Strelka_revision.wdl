version 1.0

task runStrelka {
    input {
        File normalbam
        File normalbamindex
        String tumorsample
        File tumorbam
        File tumorbamindex
        Array[File] referenceFastaFiles
    }

    File referencefasta = referenceFastaFiles[0]

    command <<<
        dir=$(echo ~{tumorbam} | sed 's/~{tumorsample}.FINAL.bam//')
        mv ~{tumorbamindex} $dir
        mv ~{normalbamindex} $dir

        /opt/strelka/bin/configureStrelkaSomaticWorkflow.py --normalBam=~{normalbam} --tumorBam=~{tumorbam} --referenceFasta=~{referencefasta} --exome --runDir=$PWD
        python2 $PWD/runWorkflow.py -m local -j 3 -g 10

        /bin/zcat < results/variants/somatic.indels.vcf.gz | awk '{if(/^##/) print; else if(/^#/) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"$0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT:"$9"\t./.:"$10"\t./.:"$11;}' - > ~{tumorsample}_Strelka_gtHeader.Somatic.snvs.vcf
        /bin/zcat < results/variants/somatic.snvs.vcf.gz | awk '{if(/^##/) print; else if(/^#/) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"$0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT:"$9"\t./.:"$10"\t./.:"$11;}' - > ~{tumorsample}_Strelka_gtHeader.Somatic.indels.vcf
    >>>

    output {
        File indelFile = "~{tumorsample}_Strelka_gtHeader.Somatic.indels.vcf"
        File snvFile = "~{tumorsample}_Strelka_gtHeader.Somatic.snvs.vcf"
    }

    runtime {
        docker: "mgibio/strelka-cwl:2.9.9"
        disks: "local-disk 400 SSD"
        memory: "32G"
        cpu: 4
    }
}