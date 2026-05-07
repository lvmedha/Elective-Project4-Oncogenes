## Print clean column headers + small head for the sheets that drive the
## GoF / pediatric-solid-tumor filter, plus full gene-level summaries
## where relevant.

suppressPackageStartupMessages({
  library(readxl)
})

xlsx <- "C:/Users/mvijayan/Downloads/ccr-24-1063_supplementary_tables_1_suppts1.xlsx"

KEY_SHEETS <- c(
  "S1 - 5009 genomic regions",
  "S7A - Overview of panels",
  "S7B - Panel gene comparsions",
  "S7C - Surrey et al Panels",
  "S7D - Reported Fusions Coverage",
  "S8 -183 drivers panel comp.",
  "S4C-Diagnostic yield SNVInDel",
  "S4D - Diagnostic yield SV",
  "S4E - Diagnostic yield CNV"
)

read_sheet_smart <- function(xlsx, sheet) {
  ## Detect title row vs header row: skip leading rows where the
  ## first cell is a long sentence (a title).
  raw <- suppressMessages(read_excel(
    xlsx, sheet = sheet, col_names = FALSE,
    col_types = "text", .name_repair = "minimal"
  ))
  first_col <- raw[[1]]
  ## Header row = first row whose first cell is short (<= 40 chars) and
  ## non-empty AND second cell is also non-empty AND not a sentence.
  hdr_row <- NA_integer_
  for (i in seq_len(min(8, nrow(raw)))) {
    v1 <- first_col[i]
    if (is.na(v1)) next
    if (nchar(v1) > 60) next
    if (ncol(raw) >= 2 && is.na(raw[[2]][i])) next
    hdr_row <- i; break
  }
  if (is.na(hdr_row)) hdr_row <- 1L
  hdr <- as.character(unlist(raw[hdr_row, ]))
  hdr[is.na(hdr) | hdr == ""] <- paste0("col", which(is.na(hdr) | hdr == ""))
  body <- raw[(hdr_row + 1L):nrow(raw), , drop = FALSE]
  colnames(body) <- hdr
  list(title_rows = if (hdr_row > 1L) raw[seq_len(hdr_row - 1L), ] else NULL,
       header = hdr,
       data = body)
}

for (s in KEY_SHEETS) {
  cat("\n\n############################################################\n")
  cat("# Sheet: ", s, "\n", sep = "")
  cat("############################################################\n")
  sm <- read_sheet_smart(xlsx, s)
  if (!is.null(sm$title_rows)) {
    cat("-- Title rows --\n")
    for (i in seq_len(nrow(sm$title_rows))) {
      vals <- as.character(unlist(sm$title_rows[i, ]))
      vals <- vals[!is.na(vals) & vals != ""]
      if (length(vals)) cat(" * ", paste(vals, collapse = " | "), "\n", sep = "")
    }
  }
  cat(sprintf("-- Header (%d cols) --\n", ncol(sm$data)))
  for (h in sm$header) cat("  - ", h, "\n", sep = "")
  cat(sprintf("-- Data dim: %d x %d --\n", nrow(sm$data), ncol(sm$data)))
  cat("-- Head (5 rows) --\n")
  show <- as.data.frame(head(sm$data, 5), stringsAsFactors = FALSE,
                        check.names = FALSE)
  show[] <- lapply(show, function(x) {
    x <- ifelse(is.na(x), "", x)
    ifelse(nchar(x) > 50, paste0(substr(x, 1, 47), "..."), x)
  })
  print(show, row.names = FALSE, right = FALSE)
}
