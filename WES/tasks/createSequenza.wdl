version 1.0

task createSequenzaFile {
    input {
        File tumorBam
        File normalBam
        File referenceFasta
        File referenceFastaIndex
        File gcWiggle
        String tumorSample
    }

    command <<<
        /usr/bin/sequenza-utils bam2seqz -n ~{normalBam} -t ~{tumorBam} --fasta ~{referenceFasta} -gc ~{gcWiggle} -o ~{tumorSample}.seqz.gz

        /usr/bin/sequenza-utils seqz_binning --seqz ~{tumorSample}.seqz.gz -w 50 -o ~{tumorSample}.binned.seqz.gz
    >>>

    output {
        File sequenzaFile = "~{tumorSample}.seqz.gz"
        File trimmedsequenza = "~{tumorSample}.binned.seqz.gz"
    }

    runtime {
        docker: "sequenza/sequenza"
        disks: "local-disk 300 SSD"
        memory: "16G"
        cpu: 1
    }

}