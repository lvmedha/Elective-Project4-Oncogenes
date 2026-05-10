11_am_missense_triage — methods (auto-generated)
==============================================

MIN_PID_RECURRENT = 2 (distinct patients at same Marker_hg38 → strong_candidate if not weak).
AM_NUM_WEAK_CUTOFF = 0.2 (if AM class missing, pathogenicity ≤ this → likely_weak).

Priority:
  1) Hotspot (Tier-0 curated codon) → strong_candidate
  2) likely_benign, or missing class with low pathogenicity → likely_weak
  3) likely_pathogenic OR recurrent (≥ MIN_PID_RECURRENT patients) → strong_candidate
  4) else → uncertain

AlphaMissense predicts damaging vs benign missense broadly; it is not a
dedicated oncogenic-activation classifier — combine with hotspots and recurrence.

Row unit: unique (gene_symbol, Marker_hg38) protein-defining variant in cohort.
