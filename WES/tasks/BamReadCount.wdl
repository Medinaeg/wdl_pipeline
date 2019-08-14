version 1.0

workflow runBRC {
    input {
        String tumor_sample
        File snv_file
        File indel_file
        File tumor_bam
        File tumor_bam_index
        File normal_bam
        File normal_bam_index
        Array[File] referenceFastaFiles
    }

    call createBED {
        input:
            tumor_sample = tumor_sample,
            snv_file = snv_file,
            indel_file = indel_file
    }

    call runBRCCommand {
        input:
            tumor_sample = tumor_sample,
            tumor_bam = tumor_bam,
            tumor_bam_index = tumor_bam_index,
            normal_bam = normal_bam,
            normal_bam_index =  normal_bam_index,
            snvBED = createBED.snvsBED,
            indelBED = createBED.indelsBED,
            referenceFastaFiles = referenceFastaFiles
    }

    output {
        Array[File] BamRCoutput = runBRCCommand.BamRCoutput
    }

}

task createBED {
    input {
        String tumor_sample
        File snv_file
        File indel_file
    }

    command <<<
        /usr/bin/vcf2bed --snvs <~{snv_file}> ~{tumor_sample}_snvs.bed

        /usr/bin/vcf2bed --deletions <~{indel_file}> ~{tumor_sample}_deletions.bed
        /usr/bin/vcf2bed --insertions <~{indel_file}> ~{tumor_sample}_insertions.bed
        
        /usr/bin/bedops --everything ~{tumor_sample}_{deletions,insertions}.bed > ~{tumor_sample}_indels.bed
    >>>

    output {
        File snvsBED = "~{tumor_sample}_snvs.bed"
        File indelsBED = "~{tumor_sample}_indels.bed"
    }

    runtime {
        docker: "biocontainers/bedops:v2.4.35dfsg-1-deb_cv1"
        disks: "local-disk 100 SSD"
        memory: "16G"
        cpu: 1
    }
}

task runBRCCommand {
    input {
        String tumor_sample
        File tumor_bam
        File tumor_bam_index
        File normal_bam 
        File normal_bam_index
        File snvBED
        File indelBED
        Array[File] referenceFastaFiles
    }

    File referenceFasta = referenceFastaFiles[0]

    command <<<
        /usr/bin/bam-readcount -q 20 -b 20 -d 500 -f ~{referenceFasta} -l ~{snvBED} -w 0 ~{tumor_bam} > ~{tumor_sample}_TUMOR_BamRC_Varscan_Strelka_snvs.txt
        /usr/bin/bam-readcount -q 20 -b 20 -d 500 -f ~{referenceFasta} -l ~{snvBED} -w 0 ~{normal_bam} > ~{tumor_sample}_NORMAL_BamRC_Varscan_Strelka_snvs.txt        

        /usr/bin/bam-readcount -q 20 -b 20 -d 500 -f ~{referenceFasta} -l ~{indelBED} -w 0 -i ~{tumor_bam} > ~{tumor_sample}_TUMOR_BamRC_Varscan_Strelka_indels.txt
        /usr/bin/bam-readcount -q 20 -b 20 -d 500 -f ~{referenceFasta} -l ~{indelBED} -w 0 -i ~{normal_bam} > ~{tumor_sample}_NORMAL_BamRC_Varscan_Strelka_indels.txt
    >>>

    output {
        Array[File] BamRCoutput = glob("*txt")
    }

    runtime {
        docker: "mgibio/bam_readcount_helper-cwl:1.1.0"
        disks: "local-disk 200 SSD"
        memory: "16G"
        cpu: 1
    }
}