## Methods: reference prep, alignment and counting

QuantSeq 3' RNA-seq reads from *Lolium rigidum* and *L. multiflorum* samples (34 libraries, single-end) were aligned to the *L. rigidum* reference genome (`Lolium_rigidum_reference.v01.fasta`, annotation `LolRi_v01.0.gff`) and counted per gene. Scripts are Slurm batch jobs (`.sb`) and are run in the order below.

### 1. Reference preparation (`prep_phytozome_inputs.sb`)

- Sequences were classified as chromosomes (seqids matching `^Chr0[0-9]+$`) or scaffolds. Scaffolds were retained (`KEEP_SCAFFOLDS=yes`) so reads mapping to them are not forced onto chromosomes. The assembly was written to `LolRi.genome.fa`.
- The GFF3 was filtered to retained sequences and cleaned with AGAT (`agat_convert_sp_gxf2gxf.pl`), which builds the missing exon features from CDS and UTR lines. Gene and mRNA counts were checked to be unchanged, and every mRNA was checked to have at least one exon.
- The cleaned annotation was converted to GTF with `agat_convert_sp_gff2gtf.pl`. Semicolons inside quoted attribute values were replaced with underscores. The result (`LolRi.star.gtf`) is used only for building the STAR index.
- For counting, each gene was extended 1,000 bp on its 3' end (`bedtools slop -s -l 0 -r 1000`, strand-aware, clipped at sequence ends) to capture reads from unannotated 3' UTR sequence typical of QuantSeq. This gave `LolRi_genes_ext1000bp.gff3` and the featureCounts annotation `LolRi_genes_ext1000bp.saf` (48,588 genes). Gene IDs were checked to agree between the SAF and the STAR GTF. The script also reports how many extensions overlap a neighbouring gene on the same or the opposite strand.

### 2. Read trimming

Raw single-end reads were trimmed with fastp (`fastp.sb`, one Slurm array task per library, 8 threads). The raw file ID is read from column 1 of `Lol_metadata.txt` and the output is named with the sample name in column 2 (`03_outputs/01_fastp/<sample>_trim.fastq.gz`). Settings:

- `--adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC` (explicit 3' adapter)
- `--trim_front1 12` (removes the first 12 bases of each read, the front-end clipping recommended for QuantSeq FWD libraries without UMIs)
- `--cut_front` and `--cut_tail` (sliding-window quality trimming from both ends, fastp default window and quality threshold)
- `--length_required 20` (reads shorter than 20 bp after trimming are discarded)
- Per-sample fastp HTML and JSON reports are written alongside the trimmed reads.

No poly-A/poly-X trimming (`--trim_poly_x`) or UMI processing/deduplication was applied.

### 3. Genome index and alignment (`star_genomegenerate.sb`, `alignment_array.sb`)

- **Index:** STAR 2.7.11b, `--runMode genomeGenerate` with `LolRi.genome.fa` and `LolRi.star.gtf`, `--sjdbOverhang 74`, `--genomeSAindexNbases 13`, 8 threads. The index is built from the unextended GTF.
- **Alignment:** One Slurm array task per sample (8 CPUs, 64 GB), with sample IDs read from column 2 of `Lol_metadata.txt`. Each trimmed single-end FASTQ was aligned with STAR 2.7.11b and these key parameters:
  - `--outSAMtype BAM SortedByCoordinate`
  - `--outFilterMultimapNmax 20`
  - `--winAnchorMultimapNmax 100`
  - `--outFilterMismatchNmax 999` with `--outFilterMismatchNoverLmax 0.08`
  - `--outFilterScoreMinOverLread 0.5` and `--outFilterMatchNminOverLread 0.5`
  - `--alignIntronMax 25000`
  - `--outFilterIntronMotifs RemoveNoncanonicalUnannotated`
  - `--outSAMattributes NH HI AS NM MD`
- BAMs were indexed with SAMtools 1.19.2. Reads aligning to more than 20 loci are discarded by STAR.

### 4. Counting (`featureCounts.sb`)

featureCounts (Subread v2.0.6) was run on all 34 BAMs, in the metadata order, against the extended SAF with the following settings. Counts are per gene (meta-feature level).

| Setting | Value |
|---|---|
| Library | single-end, stranded (`-s 1`) |
| Annotation | `-F SAF`, `LolRi_genes_ext1000bp.saf` |
| Overlaps | `-O --largestOverlap` (reads overlapping several genes are counted for the gene with the largest overlap) |
| Threads | 8 |

Two count matrices were produced:

- **`counts_unique.txt`**: uniquely mapped reads only (multi-mappers not counted). This is the primary matrix for DESeq2.
- **`counts_multi.txt`**: adds `-M --fraction`, so multi-mapping reads are split fractionally across all loci. This is used as a sensitivity check and to flag multi-copy genes.

featureCounts reports alignments rather than reads, so multi-mapping reads (up to 20 alignments each) inflate the "total alignments" and the `Unassigned_MultiMapping` category in `*.summary`. For read-level mapping rates use STAR's `*_Log.final.out`.

### 5. Downstream

Counts from `counts_unique.txt` go to DESeq2. Fractional counts are rounded before `DESeqDataSetFromMatrix` if used.

### Software

STAR 2.7.11b, SAMtools 1.19.2, Subread (featureCounts v2.0.6), AGAT 1.4.3, SeqKit 2.3.0, BEDTools 2.31.0, fastp (conda environment `fastp`; 1.3.3).