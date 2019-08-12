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

    call StrelkaCommand {
        input:
            normalbam = normalbam,
            normalbamindex = normalbamindex,
            tumorsample = tumorsample,
            tumorbam = tumorbam,
            tumorbamindex = tumorbamindex,
            referenceFastaFiles = referenceFastaFiles
    }

    call concatVCFs {
        input:
            snvFile = StrelkaCommand.snvFile,
            indelFile = StrelkaCommand.indelFile,
            tumorsample = tumorsample
    }

    output {
        File outputStrelkasnvs = concatVCFs.outputStrelkasnvs
        File outputStrelkaindels = concatVCFs.outputStrelkaindels
        File finalStrelka = concatVCFs.finalStrelka
        File finalStrelkaindex = concatVCFs.finalStrelkaIndex
    }
}

task StrelkaCommand {
    input {
        File normalbam
        File normalbamindex
        String tumorsample
        File tumorbam
        File tumorbamindex
        Array[File] referenceFastaFiles
    }

    File referenceFasta = referenceFastaFiles[0]

    command <<<
        dir=$(echo ~{tumorbam} | sed 's/~{tumorsample}.FINAL.bam//')
        mv ~{tumorbamindex} $dir
        mv ~{normalbamindex} $dir

        /opt/strelka/bin/configureStrelkaSomaticWorkflow.py --normalBam=~{normalbam} --tumorBam=~{tumorbam} --referenceFasta=~{referenceFasta} --exome --runDir=$PWD
        python2 $PWD/runWorkflow.py -m local -j 3 -g 10
        mv results/variants/somatic.indels.vcf.gz results/variants/~{tumorsample}_Strelka.somatic.indels.vcf.gz
        mv results/variants/somatic.snvs.vcf.gz results/variants/~{tumorsample}_Strelka.somatic.snvs.vcf.gz        
    >>>

    output {
        File indelFile = "results/variants/~{tumorsample}_Strelka.somatic.indels.vcf.gz"
        File snvFile = "results/variants/~{tumorsample}_Strelka.somatic.snvs.vcf.gz"
    }

    runtime {
        docker: "mgibio/strelka-cwl:2.9.9"
        disks: "local-disk 400 SSD"
        memory: "32G"
        cpu: 4
    }
}

task concatVCFs {
    input {
        File snvFile
        File indelFile
        String tumorsample
    }

    command <<<
        /usr/bin/zcat < ~{snvFile} | awk '{if(/^##/) print; else if(/^#/) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"$0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT:"$9"\t./.:"$10"\t./.:"$11;}' - > ~{tumorsample}_Strelka_gtHeader.Somatic.snvs.vcf
        /usr/bin/zcat < ~{indelFile} | awk '{if(/^##/) print; else if(/^#/) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"$0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tGT:"$9"\t./.:"$10"\t./.:"$11;}' - > ~{tumorsample}_Strelka_gtHeader.Somatic.indels.vcf

        find . -name "*.vcf" -exec /usr/local/bin/bgzip -f {} \;
        find . -name "*.vcf.gz" -exec /usr/local/bin/tabix -f {} \;

        /usr/local/bin/bcftools concat -a -o ~{tumorsample}_Strelka_gtHeader_Concat.Somatic.vcf.gz -O z ~{tumorsample}_Strelka_gtHeader.Somatic.snvs.vcf.gz ~{tumorsample}_Strelka_gtHeader.Somatic.indels.vcf.gz
        /usr/local/bin/tabix -f ~{tumorsample}_Strelka_gtHeader_Concat.Somatic.vcf.gz
    >>>

    output {
        File outputStrelkasnvs = "~{tumorsample}_Strelka_gtHeader.Somatic.snvs.vcf.gz"
        File outputStrelkaindels = "~{tumorsample}_Strelka_gtHeader.Somatic.indels.vcf.gz"
        File finalStrelka = "~{tumorsample}_Strelka_gtHeader_Concat.Somatic.vcf.gz"
        File finalStrelkaIndex = "~{tumorsample}_Strelka_gtHeader_Concat.Somatic.vcf.gz.tbi"
    }

    runtime {
        docker: "dockerbiotools/bcftools"
        disks: "local-disk 100 SSD"
        memory: "16G"
        cpu: 1
    }
}