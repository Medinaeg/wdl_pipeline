if [ "${STRANDNESS}" == "NA" ]; then
   echo hisat2 -x ${HISAT} --rg-id ${id} --rg PL:ILLUMINA --rg PU:${sample} --rg LB:${id}.${sm} --rg SM:${sample} ${files} -S ${OUT}/alignment/${sample}/align.sam
   hisat2 -x ${HISAT} --rg-id ${id} --rg PL:ILLUMINA --rg PU:${sample} --rg LB:${id}.${sm} --rg SM:${sample} ${files} -S ${OUT}/alignment/${sample}/align.sam
else
   echo hisat2 -x ${HISAT} --rg-id ${id} --rg PL:ILLUMINA --rg PU:${sample} --rg LB:${id}.${sm} --rg SM:${sample} --rna-strandness ${STRANDNESS} ${files} -S ${OUT}/alignment/${sample}/align.sam
   hisat2 -x ${HISAT} --rg-id ${id} --rg PL:ILLUMINA --rg PU:${sample} --rg LB:${id}.${sm} --rg SM:${sample} --rna-strandness ${STRANDNESS} ${files} -S ${OUT}/alignment/${sample}/align.sam
fi
