version 1.0

import "./tasks/varcalling_Mutect2.wdl" as Mutect
import "./tasks/concatVCFs.wdl" as concat

#import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/varcalling_Mutect2.wdl" as Mutect
#import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/varcalling_Varscan2.wd" as Varscan
#import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/varcalling_Strelka.wdl" as Strelka
#import "https://raw.githubusercontent.com/kcampbel/wdl_pipeline/master/tasks/varcalling_SomaticSniper.wdl" as SomaticSniper

workflow SomaticVaraintDetection{
    input {
        File fofn_bams
        File interval_list
        File reference_fasta
        File gnomad_vcf
    }

    Array[Array[String]] map_bams = read_tsv(fofn_bams)

    scatter (pair in map_bams) {
        String tumor_sample = pair[0]
        String normal_sample = pair[1]
        File tumor_bam = pair[2]
        File normal_bam = pair[3]

        call Mutect.runMutect2 as Mutect {
            input:
                interval_list = interval_list,
                reference_fasta = reference_fasta,
                gnomad_vcf = gnomad_vcf,
                tumor_sample = tumor_sample,    
                normal_sample = normal_sample,
                tumor_bam = tumor_bam,
                normal_bam = normal_bam
            }

        call concat.concatVCFs as concat {
            input:
                tumor_sample = tumor_sample,
                vcfFiles = Mutect.vcfFile
        }

    }

    output {
        Array[File] outputvcfFiles = concat.concatVCFfiles
    }
    
}

#call Mutect (tasks/varcalling_Mutect2.wdl)
# # call Varscan (tasks/varcalling_Varscan2.wdl)
# # Strelka
# # SomaticSniper
