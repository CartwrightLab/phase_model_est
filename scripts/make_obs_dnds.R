options(tidyverse.quiet = TRUE, readr.show_col_types = FALSE)
library(jsonlite)
library(tidyverse)

CSV_DIR <- "results/sampled_aln_dnds"

files <- fs::dir_ls(CSV_DIR, glob="*.csv")

tables <- map(files, read_csv, show_col_types = FALSE)

tables <- list_rbind(tables)

cat(format_csv(tables))
