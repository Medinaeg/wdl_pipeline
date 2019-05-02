version 1.0

struct ReferenceFasta {
  File reference_fasta
  File reference_fasta_index
  File reference_gtf

  String hisat_index
  File kallisto_index
}
