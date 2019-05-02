version 1.0

## This script will check to make sure that the files exist.

workflow FindFASTQs {
  input {
    String sampleName
    File fastq1
    File fastq2
  }

  command {
    if [ -f "$fastq"]
  }

  output {
    String out = read_string(stdout())
  }
}
