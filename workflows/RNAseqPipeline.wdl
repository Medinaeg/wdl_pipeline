## This WDL pipeline implements data pre-processing (trimming by Flexbar), alignment
## (by HiSat2), and expression estimation (by HTSeq for counts, Stringtie for FPKM/TPM)
## and expression estimation by pseudoalignment (by Kallisto and Pizzly). Finally,
## all expression values are aggregated into matrices for downstream analysis.
##
## Requirements/expectations :
## - One or more read groups, all belonging to a single sample
## - FASTQ files must be paired using Illumina naming properties (https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm)

# Local import
import "../tasks/FindFASTQs.wdl" as FindFiles
import "../tasks/AlignmentHisat2.wdl" as ToBam
import "../tasks/MergeSampleBams.wdl" as MergedBam
import "../tasks/ExpressionStringtie.wdl", as StringtieExpression
import "../tasks/ExpressionHTseq.wdl" as HTseqExpression
import "../tasks/ExpressionKallisto.wdl" as KallistoQuant
import "../tasks/FusionPizzly.wdl" as PizzlyFusion
import "../structs/RNAseqStructs.wdl"


# WORKFLOW DEFINITION
workflow rnaseq_workflow {
  input {
    ReferenceFasta references

    String data_directory
    String output_directory

    String? strandness
   }

   call FindFiles.GetFASTQs {

   }

   call ToBam.AlignmentHisat2 {
    input:
      references = references,


   }

   # Outputs that will be retained when execution is complete
   output {
      Array[File?] quality_yield_metrics = CollectQualityYieldMetrics.quality_yield_metrics
      File? read_group_alignment_summary_metrics = CollectReadgroupBamQualityMetrics.alignment_summary_metrics

      File? agg_alignment_summary_metrics = CollectAggregationMetrics.alignment_summary_metrics

      File output_cram = ConvertToCram.output_cram
      File output_cram_index = ConvertToCram.output_cram_index
      File output_cram_md5 = ConvertToCram.output_cram_md5

      File validate_cram_file_report = ValidateCram.report

      File output_stringtie_expression = StringtieExpression.output_stringtie_expression
      File output_kallisto_quant = KallistoQuant.output_kallisto_quant
      File output_pizzly_fusion = PizzlyFusion.output_pizzly_fusion
   }
}
