version 1.0

workflow runStrelka {
    input {
        File normalbam
        File normalbamindex
        String tumorsample
        File tumorbam
        File tumorbamindex
        Array[File] referenceFastaFiles
    }

    call runStrelkaCommand {
        input:
            normalbam = normalbam,
            normalbamindex = normalbamindex,
            tumorsample = tumorsample,
            tumorbam = tumorbam,
            tumorbamindex = tumorbamindex,
            referenceFastaFiles = referenceFastaFiles
    }

    call indexVCF {
        input:
            snvFile = runStrelkaCommand.snvFile,
            indelFile = runStrelkaCommand.indelFile,
            tumorsample = tumorsample
    }

    output {
        File indelFile = runStrelkaCommand.indelFile
        File indelFileIndex = indexVCF.indelStrelkaIndex
        File snvFile = runStrelkaCommand.snvFile
        File snvFileIndex = indexVCF.snvStrelkaIndex
    }
}

task runStrelkaCommand {
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
        File indelFile = "results/variants/~{tumorsample}_Strelka_gtHeader.Somatic.indels.vcf"
        File snvFile = "results/variants/~{tumorsample}_Strelka_gtHeader.Somatic.snvs.vcf"
    }

    runtime {
        docker: "mgibio/strelka-cwl:2.9.9"
        disks: "local-disk 400 SSD"
        memory: "32G"
        cpu: 4
    }
}

task indexVCF {
    input {
        File snvFile
        File indelFile
        String tumorsample
    }

    command <<<
        /usr/local/bin/tabix -f ~{snvFile}
        /usr/local/bin/tabix -f ~{indelFile}
    >>>

    output {
        File snvStrelkaIndex = "~{tumorsample}_Strelka_gtHeader.Somatic.snvs.vcf.tbi"
        File indelStrelkaIndex = "~{tumorsample}_Strelka_gtHeader.Somatic.indels.vcf.tbi"
    }

    runtime {
        docker: "dockerbiotools/bcftools"
        disks: "local-disk 100 SSD"
        memory: "16G"
        cpu: 1
    }
}