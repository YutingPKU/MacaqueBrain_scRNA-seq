# Generate Cellranger reference for Macaque

1. Prepare the input files
  -  Fasta files: primary sequence files from [ensembl](ftp://ftp.ensembl.org/pub/release-101/fasta/macaca_mulatta/dna/)
  `Macaca_mulatta.Mmul_10.dna.primary_assembly.Y.fa`
  - GTF files: gene annotation file from [ensembl](ftp://ftp.ensembl.org/pub/release-101/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.101.chr.gtf.gz)
  `Macaca_mulatta.Mmul_10.101.chr.gtf`

2. Filter the GTF file

  To avoid reads mapping to multiple genes, it's better to filter the gtf file.
  `cellranger mkgtf Macaca_mulatta.Mmul_10.101.chr.gtf Macaca_mulatta.Mmul_10.101.chr.filtered.gtf --attribute=gene_biotype:protein_coding --attribute=gene_biotype:lincRNA --attribute=gene_biotype:antisense`

  For single-nuclei RNA-seq, we need to count reads align to intronic regions. It's necessary to generate a 'pre-mRNA' reference.
  >listing each gene transcript locus as an exon. Thus these intronic reads will be included in the UMI counts for each gene and barcode.

  >Extract GTF annotation rows for transcripts based on the feature type transcript (column 3) of the original tab-delimited GTF and replace the feature type from transcript to exon. Here's a script to do this using the Linux utility awk.

  `awk 'BEGIN{FS="\t"; OFS="\t"} $3 == "transcript"{ print; $3="exon"; $9 = gensub("(transcript_id\\s\"{0,1})([^;\"]+)(\"{0,1});", "\\1\\2_premrna\\3;", "g", $9); print; next}{print}' Macaca_mulatta.Mmul_10.101.chr.filtered.gtf > Macaca_mulatta.Mmul_10.101.chr.filtered.premrna.gtf
`

3. Make reference package

`cellranger mkref --genome=Mmul_10.101_premrna --fasta=fasta/Macaca_mulatta.Mmul_10.dna.primary_assembly.fa --genes=Macaca_mulatta.Mmul_10.101.chr.filtered.premrna.gtf --nthreads 20 --memgb 32 --ref-version Mmul_10.101_premrna`
