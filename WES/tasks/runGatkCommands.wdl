version 1.0

workflow GatkCommands {
    
    input {
        String sample
        File duplicatesBam
        File reference_fa
        File thousG
        File knownIndels
    }
}

task BaseRecalibrate {
    input {
        String sample
        File duplicatesBam
        File reference_fa
        File thousG
        File knownIndels
        File dbsnp
    }

    command <<<
      gatk BaseRecalibrator -I ~{duplicatesBam} -R ./{reference_fa} --known-sites ~{thousG} --known-sites ~{knownIndels} --known-sites ~{dbsnp} -O ~{sample}.bqsr.table
    >>>

    output {
        File bqsrTable = "~{sample}.bqsr.table"
    }

    # runtime {

    # }
}

task ApplyBQSR {
    input {
        String sample
        File duplicatesBam
        File reference_fa
        File bqsrTable
    }

    command <<<
    gatk ApplyBQSR -R ~{reference_fa} -I ~{duplicatesBam} -bqsr ~{bqsrTable} --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 -O ~{sample}.FINAL.bam
    >>>

    output {
        File bamFinal = "~{sample}.FINAL.bam"
    }

    # runtime {

    # }

}