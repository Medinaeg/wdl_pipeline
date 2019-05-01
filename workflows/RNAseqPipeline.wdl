## This WDL pipeline implements data pre-processing (trimming by Flexbar), alignment
## (by HiSat2), and expression estimation (by HTSeq for counts, Stringtie for FPKM/TPM)
## and expression estimation by pseudoalignment (by Kallisto and Pizzly). Finally,
## all expression values are aggregated into matrices for downstream analysis.
##
## Requirements/expectations :
## - One or more read groups, all belonging to a single sample
## - FASTQ files must be paired using Illumina naming properties (https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm)

# WORKFLOW DEFINITIONN
workflow rnaseq_workflow {

   String data_directory
   String output_directory

   File adapters
   String adapter_trim_end
   Int adapter_min_overlap
   Int max_uncalled

   String hisat_index
   File reference_gtf

   File kallisto_index
   

   # Outputs that will be retained when execution is complete
   output {
      Array[File?] quality_yield_metrics = CollectQualityYieldMetrics.quality_yield_metrics
      File? read_group_alignment_summary_metrics = CollectReadgroupBamQualityMetrics.alignment_summary_metrics

      File? agg_alignment_summary_metrics = CollectAggregationMetrics.alignment_summary_metrics

      File output_cram = ConnvertToCram.output_cram
      File output_cram_index = ConvertToCram.output_cram_index
      File output_cram_md5 = ConnvertToCram.output_cram_md5

      File validate_cram_file_report = ValidateCram.report

      File output_stringtie_expression = StringtieExpression.output_stringtie_expression
      File output_kallisto_quant = KallistoQuant.output_kallisto_quant
      File output_pizzly_fusion = PizzlyFusion.output_pizzly_fusion
   }
}

if [ "${STRANDNESS}" == "NA" ]; then
   echo hisat2 -x ${HISAT} --rg-id ${id} --rg PL:ILLUMINA --rg PU:${sample} --rg LB:${id}.${sm} --rg SM:${sample} ${files} -S ${OUT}/alignment/${sample}/align.sam
   hisat2 -x ${HISAT} --rg-id ${id} --rg PL:ILLUMINA --rg PU:${sample} --rg LB:${id}.${sm} --rg SM:${sample} ${files} -S ${OUT}/alignment/${sample}/align.sam
else
   echo hisat2 -x ${HISAT} --rg-id ${id} --rg PL:ILLUMINA --rg PU:${sample} --rg LB:${id}.${sm} --rg SM:${sample} --rna-strandness ${STRANDNESS} ${files} -S ${OUT}/alignment/${sample}/align.sam
   hisat2 -x ${HISAT} --rg-id ${id} --rg PL:ILLUMINA --rg PU:${sample} --rg LB:${id}.${sm} --rg SM:${sample} --rna-strandness ${STRANDNESS} ${files} -S ${OUT}/alignment/${sample}/align.sam
fi
echo samtools view -Sb -o ${OUT}/alignment/${sample}/align.bam ${OUT}/alignment/${sample}/align.sam
samtools view -Sb -o ${OUT}/alignment/${sample}/align.bam ${OUT}/alignment/${sample}/align.sam
# All intermediate files stay in a subdirectory until the final file.
echo samtools sort -o ${OUT}/alignment/${sample}_FINAL.bam ${OUT}/alignment/${sample}/align.bam
samtools sort -o ${OUT}/alignment/${sample}_FINAL.bam ${OUT}/alignment/${sample}/align.bam
rm ${OUT}/alignment/${sample}/align.bam
/u/flashscratch/k/katiecam/software/stringtie-1.3.5/stringtie -G ${GTF} -e -B -o ${OUT}/results/${sample}_transcripts.gtf -A ${OUT}/results/${sample}_abundances.tsv ${OUT}/alignment/${sample}_FINAL.bam

/u/flashscratch/k/katiecam/miniconda2/bin/htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id ${OUT}/alignment/${sample}_FINAL.bam ${GTF} > ${OUT}/results/${sample_FINAL}_HTseqCounts.tsv

/u/flashscratch/k/katiecam/miniconda2/bin/kallisto quant -i ${KALLISTO} -b 100 --fusion -o ${OUT}/results/${sample} ${files}

/u/flashscratch/k/katiecam/miniconda2/bin/pizzly -k 31 --gtf ${GTF} --align-score 2 --insert-size 400 --fasta ${CDNA} --output ${OUT}/results/${sample}/pizzly ${OUT}/results/${sample}/fusion.txt
