# Snappy: *de novo* identifiaction of methylation sites

Snappy provides fast and accurate identification of methylation sites based on 
Oxford Nanopore reads. Snappy combines a new-graph based enrichment 
algorithm with simultaneous analysis of the raw sequencing data, that 
significantly enhances accuracy of motifs identification.
Snappy is mainly designed to work with data obtained from ONT r10.4.1 flow cells.

## Installation

Snappy is availabe in PyPi:

```
(base) $ conda create -n snappy python=3.12
(base) $ conda activate snappy
(snappy) $ pip install ont-snappy
```

## Usage

Let's suppose that we have pod5/ folder with raw whole genome sequencing data of the target OBJECT.
Typical data processing pipeline including Snappy will look as follows:
1. Rebasecalling of the raw POD5 files using special modified basecalling models:
```
dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v5.0.0 pod5/ --modified-bases-models dna_r10.4.1_e8.2_400bps_sup@v5.0.0_6mA@v1,dna_r10.4.1_e8.2_400bps_sup@v5.0.0_4mC_5mC@v1 > OBJECT.bam
```
2. Obtaining genome assembly with Flye, if it's nessesary:
```
samtools fastq OBJECT.bam | flye -t 32 -o OBJECT_assembly -nano-raw -
```

3. Mapping the reads to the genome assembly:
```
samtools fastq -T MM,ML OBJECT.bam | minimap2 -ax map-ont -y OBJECT_assembly/assembly.fasta -t 32 - | samtools view -b | samtools sort > OBJECT.mapped.bam
samtools index OBJECT.mapped.bam
```

4. Running Modkit:
```
modkit pileup -t 32 OBJECT.mapped.bam OBJECT_modkit.bed -r OBJECT_assembly --only-tabs --filter-threshold 0.66
```

5. Running Snappy:
```
snappy -mk_bed OBJECT_modkit.bed -genome OBJECT_assembly/assembly.fasta -outdir OBJECT_snappy
```


## Snappy output

The main Snappy output is two text files `Summury_table.txt` and `Results_table.txt`. The first file presents summarized information about all identified motifs. The second file provides methylation profiles for each identified motif divided by contigs. For more advanced users, Snappy saves extended and filtered modkit table used for the enrichment, and regexp-formatted records for each identified motif so to provide convinient access to the data using Polars or Pandas. Finally, the output directory `VIZ` contians vizualization for each identified motif.

## For MicrobeMod users

If you have already processed your data with MicrobeMod, you can run Snappy using *_low.bed file provided by MicrobeMod. For example:

```
snappy -mk_bed [MicrobeMod output]/*_low.bed -genome [genome.fasta] -outdir snappy_out
```

## Citation

Dmitry N Konanov, Danil V Krivonos, Validislav V Babenko, Elena N Ilina. **Snappy: *de novo* identification of DNA methylation sites based on Oxford Nanopore reads**. *bioRxiv* 2025.08.03.668330; [Link to preprint](https://doi.org/10.1101/2025.08.03.668330)