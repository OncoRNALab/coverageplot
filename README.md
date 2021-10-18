# coverageplot

Welcome! This script generates an interactive figure of the coverage of a specified genomic region in a way similar to IGV. You can start from a BAM file or a BED file. Additionally, you can submit LNA or ASO sequences to visualise where these oligonucleotides would bind based on sequence compelementarity. Please let me know if there are anything else you would like added to the plots or if you have other suggestions. 

## Dependencies
For coverageplot to work, it is important you have the following software installed:
  - Bedtools
  - Samtools
  - R and R packages: tidyverse, plotly, ggpubr, htmlwidgets, spgs and optparse

## Use
You can get more information of the different parameters by running the script with -h: 
`Rscript coverageplot.R -h`

This will display: 
```
Usage: coverageplot.R [options]
Options:
	-r CHARACTER, --region=CHARACTER
		Region of the genome that you want to visualise. Can either be a region copied from Ensembl or the Ensembl gene id

	-b CHARACTER, --bedbam=CHARACTER
		BED or BAM file with the correct extension

	-g CHARACTER, --gtf=CHARACTER
		GTF file

	-f CHARACTER, --fasta=CHARACTER
		Whole genome FASTA file

	-l CHARACTER, --lna=CHARACTER
		If you want to visualise where an LNA or ASO will bind, use this argument. If multiple sequences, separate by ','.

	-h, --help
		Show this help message and exit
```

The output of the script is an interactive html file showing the region of interest. 

## Example
An example of a command to start the program would be: 
```
Rscript --vanilla coverageplot.R -r "Chromosome MT: 1,671-3,229" -f genome.fa -g Homo_sapiens.GRCh38.91_withspikes.gtf -b RNA016611_S31_R1_001_Aligned.sortedByCoord.out.bam -l "AAGAGCACACCCGTCT,CCAACACAGGCATGCT"
```

This command will run the script for the genomic region of MT-RNR2 (note: I could also have used ENSG00000210082) starting from a BAM file (note: you can also start from a bed file, just make sure the files have the right extensions. This specific command will also display ASO sequences AAGAGCACACCCGTCT and CCAACACAGGCATGCT. The script allows no errors in the sequences and will also look for the reverse complement. 
