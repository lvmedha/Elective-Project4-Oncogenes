08_vep_input_full.* -- VEP-ready inputs (whole VAF)
===================================================

Source : C:/Users/mvijayan/Documents/AllMarkers_VAF_long.tsv
Variants: 4251 unique (chrom, pos, ref, alt) on GRCh38
Genome  : GRCh38 / hg38, UCSC chromosome naming (chr1..chr22, chrX, chrY)
Scope   : whole cohort, panel-agnostic. AlphaMissense will only score
          missense substitutions; non-missense rows will get a normal VEP
          annotation but an empty AM score, which is expected.

Files (this folder):
  08_vep_input_full.vcf       VCFv4.2 (use with `vep -i ... .vcf`)
  08_vep_input_full.ensembl   VEP default tab format (Part 3.1):
                              chrom start end allele strand identifier
                              (use with `vep -i ... .ensembl --format ensembl`)

Joining AM scores back into the project
---------------------------------------
Once an annotated VCF is produced (vep + --plugin AlphaMissense), parse
the CSQ INFO field, emit a TSV keyed by (chrom, pos, ref, alt) carrying
am_pathogenicity / am_class / LoF / gene / consequence, and join onto
results/ped_gof_snv/03_vaf_annotated.tsv (or results/<RUN_NAME>/ for another panel) by
(chrom, pos, ref, alt) -- same key as this VCF.

Notes
-----
  - chrom naming: this file uses chr1..chr22, chrX, chrY (UCSC-style)
    to match a UCSC-style VEP cache. For an Ensembl-style cache,
    strip the chr prefix or pass --synonyms to vep (chr_synonyms.txt path).
  - For variants where the AM TSV is missing a row (synonymous, intronic,
    indels longer than 1 nt), VEP writes the row with am_pathogenicity
    blank -- expected, not an error.
