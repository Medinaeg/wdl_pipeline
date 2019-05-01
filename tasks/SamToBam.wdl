samtools view -Sb -o ${OUT}/alignment/${sample}/align.bam ${OUT}/alignment/${sample}/align.sam
samtools sort -o ${OUT}/alignment/${sample}_FINAL.bam ${OUT}/alignment/${sample}/align.bam
