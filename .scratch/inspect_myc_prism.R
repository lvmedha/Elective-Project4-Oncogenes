x <- data.table::fread("results/ped_gof_snv/10_prism_mut_vs_wt_drug_tests.tsv")
m <- x[mutated_gene == "MYC"]
cat("MYC PRISM tests\n")
cat("  total drugs        :", nrow(m), "\n")
cat("  n_mut range        :", paste(range(m$n_mut), collapse = "-"), "\n")
cat("  n_wt  range        :", paste(range(m$n_wt),  collapse = "-"), "\n")
cat("  min p-value        :", signif(min(m$p_wilcox, na.rm = TRUE), 3), "\n")
cat("  min FDR            :", signif(min(m$fdr_within_gene, na.rm = TRUE), 3), "\n")
cat("  drugs FDR<0.05     :", sum(m$fdr_within_gene < 0.05, na.rm = TRUE), "\n")
cat("  drugs p<0.001      :", sum(m$p_wilcox < 0.001, na.rm = TRUE), "\n")
cat("  drugs p<0.01       :", sum(m$p_wilcox < 0.01, na.rm = TRUE), "\n\n")
cat("Top 8 MYC drugs by p-value:\n")
print(m[order(p_wilcox)][1:8, .(CompoundName, TargetOrMechanism,
                                delta = round(delta_mut_minus_wt, 3),
                                p = signif(p_wilcox, 3),
                                fdr = signif(fdr_within_gene, 3))])
