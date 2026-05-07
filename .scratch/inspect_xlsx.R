## Inspect ccr-24-1063 supplementary tables: list sheets, dims, and head rows.
suppressPackageStartupMessages({
  if (!requireNamespace("readxl", quietly = TRUE)) {
    install.packages("readxl", repos = "https://cloud.r-project.org")
  }
  library(readxl)
})

xlsx <- "C:/Users/mvijayan/Downloads/ccr-24-1063_supplementary_tables_1_suppts1.xlsx"
stopifnot(file.exists(xlsx))

sheets <- excel_sheets(xlsx)
cat("=== Sheets (", length(sheets), ") ===\n", sep = "")
for (s in sheets) cat(" - ", s, "\n", sep = "")

for (s in sheets) {
  cat("\n\n========================================\n")
  cat("Sheet: ", s, "\n", sep = "")
  cat("========================================\n")
  df <- tryCatch(
    suppressMessages(read_excel(xlsx, sheet = s, col_names = FALSE,
                                col_types = "text", .name_repair = "minimal")),
    error = function(e) { cat("READ ERROR:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(df)) next
  cat(sprintf("dim: %d rows x %d cols\n", nrow(df), ncol(df)))
  nshow <- min(20L, nrow(df))
  cap <- function(x, n = 80) {
    x <- ifelse(is.na(x), "", x)
    ifelse(nchar(x) > n, paste0(substr(x, 1, n - 3), "..."), x)
  }
  show <- as.data.frame(lapply(df[seq_len(nshow), , drop = FALSE], cap),
                        stringsAsFactors = FALSE,
                        check.names = FALSE)
  print(show, row.names = FALSE, right = FALSE)
}
