## This WDL pipeline implements data pre-processing (trimming by Flexbar), alignment
## (by HiSat2), and expression estimation (by HTSeq for counts, Stringtie for FPKM/TPM)
## and expression estimation by pseudoalignment (by Kallisto and Pizzly). Finally,
## all expression values are aggregated into matrices for downstream analysis.
##
## Requirements/expectations :
## - 3-column tab-delimited file, with the following format:
## -- Column 1. Sample name
## -- Column 2. Read 1 FASTQ file
## -- Column 3. Read 2 FASTQ file (Leave blank if reads are unpaired)
##

# Local import
import "../tasks/AlignmentHisat2.wdl" as FastqToAlignedBam
import "../structs/RNAseqStructs.wdl"


# WORKFLOW DEFINITION
workflow rnaseq_workflow {
  input {
    File fofn # file of file names

    ReferenceFasta references

    String? strandness
   }

   # TODO: Implement QC - Find files
   # Read in FOFN to get sample:fastq mapping
   Array[Array[String]] inputSamples = read_tsv(fofn)

   # Run alignments on all samples
   scatter (sample in inputSamples) {
    call FastqToAlignedBam.Hisat2Alignment {
      input:
        sampleName = inputSamples[0],
        fastq1 = inputSamples[1],
        fastq2 = inputSamples[2],
        strandness = strandness
    }
   }

   # Outputs that will be retained when execution is complete
   output {

   }
}
