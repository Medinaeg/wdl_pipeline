# Varscan2: https://github.com/genome/analysis-workflows/blob/master/definitions/tools/varscan_somatic.cwl
version 1.0

workflow runVarscan {
    input {
        File tumor_bam
        File tumor_bam_index
        String tumor_sample
        File normal_bam
        File normal_bam_index
        Array[File] referenceFastaFiles
    }

    call VarscanCommand {
        input:
            tumor_bam = tumor_bam,
            tumor_bam_index = tumor_bam_index,
            tumor_sample = tumor_sample,
            normal_bam = normal_bam,
            normal_bam_index = normal_bam_index,
            referenceFastaFiles = referenceFastaFiles
    }

    call compressVCFs {
        input:
        tumor_sample = tumor_sample,
        Varscanoutput = VarscanCommand.Varscanoutput

    }
    
    call VariantFilter {
        input:
        referenceFastaFiles = referenceFastaFiles,
        snpFile = compressVCFs.snpCompressed,
        snphcFile = compressVCFs.snphcCompressed,
        indelFile = compressVCFs.indelCompressed,
        indelhcFile = compressVCFs.indelhcCompressed, 
        snpIndex = compressVCFs.snpIndex,
        snphcIndex = compressVCFs.snphcIndex,
        indelIndex = compressVCFs.indelIndex, 
        indelhcIndex = compressVCFs.indelhcIndex, 
        tumor_sample = tumor_sample
    }

    call concatVCFs {
        input:
        FilterSnpFile = VariantFilter.FilterSnpFile,
        FilterSnpFileIndex = VariantFilter.FilterSnpFileIndex,
        FilterIndelFile = VariantFilter.FilterIndelFile,
        FilterIndelFileIndex = VariantFilter.FilterIndelFileIndex,
        tumor_sample = tumor_sample
    }

    output {
        Array[File] allVarscanFiles = VarscanCommand.Varscanoutput
        File finalVarscan = concatVCFs.finalVarscan
        File finalVarscanindex = concatVCFs.finalVarscanindex
    }
}

task VarscanCommand {
    input {
        File tumor_bam
        File tumor_bam_index
        String tumor_sample
        File normal_bam
        File normal_bam_index
        Array[File] referenceFastaFiles
    }

    File reference_fasta = referenceFastaFiles[0]

    command <<<
    java -jar /opt/varscan/VarScan.jar somatic \
    <(/opt/samtools/bin/samtools mpileup --no-baq -f ~{reference_fasta} ~{normal_bam} ~{tumor_bam}) \
    ~{tumor_sample}_main_output_Varscan \
    --strand-filter 0 \
    --min-coverage 20 \
    --min-var-freq 0.05 \
    --p-value 0.99 \
    --mpileup 1 \
    --output-vcf 1

    java -jar /opt/varscan/VarScan.jar processSomatic ~{tumor_sample}_main_output_Varscan.snp.vcf
    java -jar /opt/varscan/VarScan.jar processSomatic ~{tumor_sample}_main_output_Varscan.indel.vcf
    
    mv ~{tumor_sample}_main_output_Varscan.indel.Germline.hc.vcf ~{tumor_sample}_Varscan.indel.Germline.hc.vcf
    mv ~{tumor_sample}_main_output_Varscan.indel.Germline.vcf ~{tumor_sample}_Varscan.indel.Germline.vcf
    mv ~{tumor_sample}_main_output_Varscan.indel.LOH.hc.vcf ~{tumor_sample}_Varscan.indel.LOH.hc.vcf
    mv ~{tumor_sample}_main_output_Varscan.indel.LOH.vcf ~{tumor_sample}_Varscan.indel.LOH.vcf
    mv ~{tumor_sample}_main_output_Varscan.indel.Somatic.hc.vcf ~{tumor_sample}_Varscan.indel.Somatic.hc.vcf
    mv ~{tumor_sample}_main_output_Varscan.indel.Somatic.vcf ~{tumor_sample}_Varscan.indel.Somatic.vcf
    mv ~{tumor_sample}_main_output_Varscan.snp.Germline.hc.vcf ~{tumor_sample}_Varscan.snp.Germline.hc.vcf
    mv ~{tumor_sample}_main_output_Varscan.snp.Germline.vcf ~{tumor_sample}_Varscan.snp.Germline.vcf
    mv ~{tumor_sample}_main_output_Varscan.snp.LOH.hc.vcf ~{tumor_sample}_Varscan.snp.LOH.hc.vcf
    mv ~{tumor_sample}_main_output_Varscan.snp.LOH.vcf ~{tumor_sample}_Varscan.snp.LOH.vcf
    mv ~{tumor_sample}_main_output_Varscan.snp.Somatic.hc.vcf ~{tumor_sample}_Varscan.snp.Somatic.hc.vcf
    mv ~{tumor_sample}_main_output_Varscan.snp.Somatic.vcf ~{tumor_sample}_Varscan.snp.Somatic.vcf
    >>>

    output {
        Array[File] Varscanoutput = glob("*vcf*")
    }

    runtime {
        docker: "mgibio/cle:v1.3.1"
        disks: "local-disk 300 SSD"
        memory: "16G"
        cpu: 1
    }
}

task compressVCFs {
    input {
        String tumor_sample
        Array[File] Varscanoutput
    }

    File snpFile = Varscanoutput[11]
    File snphcFile = Varscanoutput[10]
    File indelFile = Varscanoutput[4]
    File indelhcFile = Varscanoutput[5]

    command <<<
        /usr/local/bin/bgzip -f ~{snpFile}
        /usr/local/bin/bgzip -f ~{snphcFile}
        /usr/local/bin/bgzip -f ~{indelFile}
        /usr/local/bin/bgzip -f ~{indelhcFile}
        dir=$(echo ~{snpFile} | sed 's/~{tumor_sample}_Varscan.snp.Somatic.vcf//')
        mv $dir/* .

        find . -name '*.vcf.gz' -exec /usr/local/bin/tabix -f {} \;
    >>>

    output {
        File snpCompressed = "~{tumor_sample}_Varscan.snp.Somatic.vcf.gz"
        File snphcCompressed = "~{tumor_sample}_Varscan.snp.Somatic.hc.vcf.gz"
        File indelCompressed = "~{tumor_sample}_Varscan.indel.Somatic.vcf.gz"
        File indelhcCompressed = "~{tumor_sample}_Varscan.indel.Somatic.hc.vcf.gz"
        File snpIndex = "~{tumor_sample}_Varscan.snp.Somatic.vcf.gz.tbi"
        File snphcIndex = "~{tumor_sample}_Varscan.snp.Somatic.hc.vcf.gz.tbi"
        File indelIndex = "~{tumor_sample}_Varscan.indel.Somatic.vcf.gz.tbi"
        File indelhcIndex = "~{tumor_sample}_Varscan.indel.Somatic.hc.vcf.gz.tbi"
    }

    runtime {
        docker: "dockerbiotools/bcftools"
        disks: "local-disk 100 SSD"
        memory: "16G"
        cpu: 1
    }
}

task VariantFilter {
    input {
        Array[File] referenceFastaFiles
        File snpFile
        File snphcFile
        File indelFile
        File indelhcFile
        File snpIndex
        File snphcIndex
        File indelIndex
        File indelhcIndex
        String tumor_sample
    }

    File referenceFasta = referenceFastaFiles[0]

    command <<<
        /usr/gitc/gatk4/gatk-launch VariantFiltration -R ~{referenceFasta} -V ~{snpFile} --mask ~{snphcFile} --maskName "processSomatic" --filterNotInMask -O ~{tumor_sample}_Varscan_Gatk_filtered.snp.Somatic.hc.vcf.gz
        /usr/gitc/gatk4/gatk-launch VariantFiltration -R ~{referenceFasta} -V ~{indelFile} --mask ~{indelhcFile} --maskName "processSomatic" --filterNotInMask -O ~{tumor_sample}_Varscan_Gatk_filtered.indel.Somatic.hc.vcf.gz
    >>>

    output {
        File FilterSnpFile = "~{tumor_sample}_Varscan_Gatk_filtered.snp.Somatic.hc.vcf.gz"
        File FilterSnpFileIndex = "~{tumor_sample}_Varscan_Gatk_filtered.snp.Somatic.hc.vcf.gz.tbi"
        File FilterIndelFile = "~{tumor_sample}_Varscan_Gatk_filtered.indel.Somatic.hc.vcf.gz"
        File FilterIndelFileIndex = "~{tumor_sample}_Varscan_Gatk_filtered.indel.Somatic.hc.vcf.gz.tbi"
    }
    
    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
        disks: "local-disk 100 SSD"
        memory: "16G"
        cpu: 1
    }
}

task concatVCFs {
    input {
        File FilterSnpFile
        File FilterSnpFileIndex
        File FilterIndelFile
        File FilterIndelFileIndex
        String tumor_sample
    }

    command <<<
        /usr/local/bin/bcftools concat -a -o ~{tumor_sample}_Varscan_Gatk_filtered_Concat.Somatic.vcf.gz -O z ~{FilterSnpFile} ~{FilterIndelFile}
        /usr/local/bin/tabix -f ~{tumor_sample}_Varscan_Gatk_filtered_Concat.Somatic.vcf.gz
    >>>

    output {
        File finalVarscan = "~{tumor_sample}_Varscan_Gatk_filtered_Concat.Somatic.vcf.gz"
        File finalVarscanindex = "~{tumor_sample}_Varscan_Gatk_filtered_Concat.Somatic.vcf.gz.tbi"
    }

    runtime {
        docker: "dockerbiotools/bcftools"
        disks: "local-disk 100 SSD"
        memory: "16G"
        cpu: 1
    }
}