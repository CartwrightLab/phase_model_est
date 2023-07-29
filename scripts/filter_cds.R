library(tidyverse)
library(seqinr)

ARGS <- commandArgs(trailingOnly = TRUE)

output_dir <- ARGS[1]
input_dir <- ARGS[2]
globs <- ARGS[c(-1,-2)]

input_dir <- fs::path_abs(input_dir)

fasta_files <- globs |> map(\(x) fs::dir_ls(path=input_dir, glob=x))
fasta_files <- list_c(fasta_files)

fs::dir_create(output_dir)

fasta_files |> walk(.progress = TRUE, function(filename) {
    dna <- read.fasta(filename, as.string=TRUE)
    # convert to upper case and strip gaps and terminal stop codons
    dna <- dna |> map_chr(function(x) {
        x <- x |> str_flatten() |> str_to_upper() |>
            str_remove_all(fixed("-")) |>
            str_remove("(TAA|TAG|TGA)$")
        x
    })
    # skip any pairs that contain "N"
    if(any(str_detect(dna, "[^ACGT]"))) {
        return()
    }
    # skip any pairs that are not multiples of 3
    if(any((nchar(dna) %% 3) != 0)) {
        return()
    }
    # skip any pairs that are too long
    if(any(nchar(dna) > 9000L)) {
        return()
    }
    # skip any pairs that have internal stop codons
    b <- dna |> strsplit(character(0L)) |> map_lgl(function(x) {
        s <- seq(1, length(x), 3)

        any(paste0(x[s], x[s+1], x[s+2]) %in% c("TAA", "TAG", "TGA"))
    })
    if(any(b)) {
        return()
    }
    
    # sort sequences randomly based on entropy of the dna sequences
    h <- rlang::hash(sort(dna))
    h <- strtoi(substr(h, 1, 7), 16)
    dna <- withr::with_seed(h, as.list(sample(dna)))

    # save sequences
    basename <- fs::path_file(filename) |> fs::path_ext_remove()
    fasta <- fs::path(output_dir, basename, ext = "fasta")
    write.fasta(dna, names(dna), file.out=fasta)
})


