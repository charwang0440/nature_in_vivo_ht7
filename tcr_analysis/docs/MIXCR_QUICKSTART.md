# MiXCR Quickstart (Human TCR, full length)

## Install & License
1. Get a license from MiXCR website.
2. Download/install MiXCR.

## Analyze
Example (Takara Human TCR V2 full-length):
```
mixcr analyze takara-human-tcr-V2-full-length \
  /path/to/m9_R1.fastq.gz \
  /path/to/m9_R2.fastq.gz \
  m9-result
```

## QC reports
```
mixcr exportQc chainUsage --hide-non-functional *.clns chainUsage.pdf
mixcr exportQc align *.vdjca alignQc.pdf
```