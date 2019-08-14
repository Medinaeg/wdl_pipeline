task runSelectVariant {
    input {
        File mergedVCF
        String tumorSample
        Array[File] referenceFastaFiles
    }

    File referenceFasta = referenceFastaFiles[0]

    command <<<
        /usr/gitc/gatk4/gatk-launch SelectVariants -R ~{referenceFasta} -V ~{mergedVCF} -select 'set=="Intersection"' -O ~{tumorSample}_Merged_Intersection_indels.vcf
    >>>

    output {
        File indelvcfintersection = "~{tumorSample}_Merged_Intersection_indels.vcf"
    }

    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        disks: "local-disk 100 SSD"
        memory: "16G"
        cpu: 1
    }
}
