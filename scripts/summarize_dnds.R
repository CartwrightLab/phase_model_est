options(tidyverse.quiet = TRUE)
library(jsonlite)
library(tidyverse)

read_fasta <- function(path) {
    fasta <- as.character(ape::read.FASTA(path))
    fasta <- map_chr(fasta, str_flatten)
    fasta[] <- str_to_upper(fasta)

    seqinr::as.alignment(length(fasta), names(fasta), fasta, com = NULL)
}

process_file <- function(path) {
    fasta <- read_fasta(path)

    tab <- seqinr::kaks(fasta, verbose = TRUE)
    tab <- as_tibble(tab)
    tab <- add_column(tab, file = fs::path_file(path), .before = 1L)

    tab
}

## MAIN ########################################################################

summarize_dnds_main <- function(dir) {
    files <- fs::dir_ls(dir, type="file", glob="*.fasta")

    dat <- map(files, process_file, .progress = TRUE) |> 
        list_rbind()

    dat <- dat |> add_column(species = dir, .before = 1L)

    dat
}

## RSCRIPT BODY ################################################################

if(!rlang::is_interactive()) {
    args <- commandArgs(trailingOnly = TRUE)
    dat <- summarize_aln_main(args[1])
    cat(format_csv(dat))
}
