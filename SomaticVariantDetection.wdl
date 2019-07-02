version 1.0

#import "./tasks/varcalling_Mutect2.wdl" as Mutect
#import "./tasks/concatVCFs.wdl" as concat
#import "./tasks/varcalling_Varscan2.wdl" as Varscan

#import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/varcalling_Mutect2.wdl" as Mutect
import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/varcalling_Varscan2.wdl" as Varscan
#import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/varcalling_Strelka.wdl" as Strelka
#import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/varcalling_SomaticSniper.wdl" as SomaticSniper

workflow SomaticVaraintDetection{
    input {
        File fofn_bams
        #File interval_list
        File reference_fasta
        File reference_fasta_index
        #File gnomad_vcf
    }

    Array[Array[String]] map_bams = read_tsv(fofn_bams)

    scatter (pair in map_bams) {
        String tumor_sample = pair[0]
        String normal_sample = pair[1]
        File tumor_bam = pair[2]
        File tumor_bam_index = pair[3]
        File normal_bam = pair[4]
        File normal_bam_index = pair[5]

        # call Mutect.runMutect2 as Mutect {
        #     input:
        #         interval_list = interval_list,
        #         reference_fasta = reference_fasta,
        #         gnomad_vcf = gnomad_vcf,
        #         tumor_sample = tumor_sample,    
        #         normal_sample = normal_sample,
        #         tumor_bam = tumor_bam,
        #         normal_bam = normal_bam
        #     }

        # call concat.concatVCFs as concat {
        #     input:
        #         tumor_sample = tumor_sample,
        #         vcfFiles = Mutect.vcfFile
        # }

        call Varscan.runVarscan2 as Varscan {
            input:
                reference_fasta = reference_fasta,
                reference_fasta_index =reference_fasta_index,
                tumor_bam = tumor_bam,
                tumor_bam_index = tumor_bam_index,
                tumor_sample = tumor_sample,
                normal_bam = normal_bam,
                normal_bam_index = normal_bam_index
        }

    }

    output {
    #    Array[File] outputvcfFiles = concat.concatVCFfiles
        Array[File] outputsnpFiles = Varscan.snpFile
        Array[File] outputindelFiles = Varscan.indelFile
    }
    
}

#call Mutect (tasks/varcalling_Mutect2.wdl)
#call Varscan (tasks/varcalling_Varscan2.wdl)
# # Strelka
# # SomaticSniper
