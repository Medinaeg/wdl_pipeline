version 1.0

workflow Mutectwf {
    input {
        String sample
        File tumorbamFile
        File normalbamFile
        String tumorFileName
        String normalFileName
        File referenceFile
        Array[File] IntervalFilesList
    }

    scatter(Interval in IntervalFilesList){
        File intervalFile = IntervalFilesList[0]
        Int j = IntervalFilesList[1]

        call Mutect2 {
            input:
                sample = sample,
                intervalFile = intervalFile,
                tumorFileName = tumorFileName,
                normalFileName = normalFileName,
                tumorbamFile = tumorbamFile,
                normalbamFile = normalbamFile,
                referenceFile = referenceFile,
                j = j
        }

        call MergeVCFs {
            input:
                    sample = sample,
                    vcfFiles = Mutect2.vcfFile
        }
    }

    output {
        Array[File] mergedVCFfiles = MergeVCFs.mergedVCFfiles
    }

}

task Mutect2 {
    input {
        String sample
        File intervalFile
        String tumorFileName
        String normalFileName
        File tumorbamFile
        File normalbamFile
        File referenceFile
        Int j
    }

    command <<<     
    /usr/gitc/gatk4/gatk-launch Mutect2 -R ~{referenceFile} -I ~{tumorbamFile} -tumor ~{tumorFileName} -I ~{normalbamFile} -normal ~{normalFileName} -O ~{sample}.~{j}.vcf.gz -L ~{intervalFile}
    >>>

    output{
        File vcfFile = "~{sample}.~{j}.vcf.gz"
    }

    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        disks: "local-disk 500 SSD"
        memory: "16G"
        cpu: 1  
    }

}

task MergeVCFs {
    input{
        String sample
        Array[File] vcfFiles
    }

    command <<<
    
    >>>

    output {
        File mergedVCFfiles = "~{sample}.FINAL.mergedVCF.vcf.gz"
    }

    runtime {
        docker: "opengenomics/vcftools-tools:latest"
        disks: "local-disk 100 SSD"
        memory: "16G"
        cpu: 1  
    }
}