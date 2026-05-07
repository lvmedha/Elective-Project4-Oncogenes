08_vep_input.vcf / 08_vep_input.tsv -- VEP-ready inputs
========================================================

Source : C:/Users/mvijayan/Documents/Elective-Project4-Oncogenes/results/all_42genes/05_mutations_tiered.tsv
Variants: 268 unique (chrom, pos, ref, alt) on GRCh38
Genome  : GRCh38 / hg38, Ensembl chromosome naming (no 'chr' prefix)

Files:
  08_vep_input.vcf  -- VCFv4.2 with INFO fields GENE / HGVSP / MUT_STATUS / VCLASS / N_PIDS
  08_vep_input.tsv  -- VEP default format (CHROM POS . REF ALT . . .)

Optional strict normalisation (recommended; run inside WSL):
  bcftools norm -f /path/to/GRCh38.fa -c w 08_vep_input.vcf \
    -Oz -o 08_vep_input.norm.vcf.gz
  tabix -p vcf 08_vep_input.norm.vcf.gz

Local VEP + AlphaMissense plugin (run inside the alpha_tools conda env):
  vep \
    --input_file  08_vep_input.vcf \
    --output_file 08_vep_alphamissense.tsv \
    --species homo_sapiens --assembly GRCh38 \
    --cache --offline --dir_cache /path/to/vep_cache \
    --fasta /path/to/GRCh38.fa \
    --tab --everything --pick \
    --plugin AlphaMissense,file=/path/to/AlphaMissense_hg38.tsv.gz

Notes:
  - The AlphaMissense TSV must be bgzipped + tabix-indexed:
      tabix -s 1 -b 2 -e 2 -f -S 1 AlphaMissense_hg38.tsv.gz
  - --pick keeps one consequence per variant; drop it for all transcripts.
  - --everything turns on canonical, mane, hgvs, protein, biotype, etc.
  - If your VEP cache uses 'chr1' style names, add `--chr_synonyms` or
    rename contigs first (bcftools annotate --rename-chrs).
