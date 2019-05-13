version 1.0

workflow MarkDuplicates {
    input {
        String sample
        File bamFile
    }

    call PicardCommand {
        input:
            sample = sample,
            bamFile = bamFile
    }

    output {
        #Array[File] duplicatesBam
        #Array[File] metricFile       
    }

}


task PicardCommand {
    input {
        String sample
        File bamFile
    }

    command <<< 
        java -Xmx1g -jar /u/local/apps/picard-tools/current/picard.jar MarkDuplicates I=~{bamFile} O=~{sample}.md.bam ASSUME_SORT_ORDER=coordinate METRICS_FILE=~{sample}.md.txt QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT
    >>>

    output {
        File duplicatesBam = "~{sample}.md.bam"
        File metricFile = "~{sample}.md.txt"
    }
    
    # runtime {

    # }
}
