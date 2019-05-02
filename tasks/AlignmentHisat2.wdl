version 1.0

import "../"
import "../structs/RNAseqStructs.wdl"

workflow Hisat2Alignment {
  input {
    String sampleName
    File fastq1
    File? fastq2
    String? strandness
  }

  call
}


if [ "${STRANDNESS}" == "NA" ]; then
   echo hisat2 -x ${HISAT} --rg-id ${id} --rg PL:ILLUMINA --rg PU:${sample} --rg LB:${id}.${sm} --rg SM:${sample} ${files} -S ${OUT}/alignment/${sample}/align.sam
   hisat2 -x ${HISAT} --rg-id ${id} --rg PL:ILLUMINA --rg PU:${sample} --rg LB:${id}.${sm} --rg SM:${sample} ${files} -S ${OUT}/alignment/${sample}/align.sam
else
   echo hisat2 -x ${HISAT} --rg-id ${id} --rg PL:ILLUMINA --rg PU:${sample} --rg LB:${id}.${sm} --rg SM:${sample} --rna-strandness ${STRANDNESS} ${files} -S ${OUT}/alignment/${sample}/align.sam
   hisat2 -x ${HISAT} --rg-id ${id} --rg PL:ILLUMINA --rg PU:${sample} --rg LB:${id}.${sm} --rg SM:${sample} --rna-strandness ${STRANDNESS} ${files} -S ${OUT}/alignment/${sample}/align.sam
fi

call SamToBam
