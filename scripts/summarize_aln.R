options(tidyverse.quiet = TRUE)
library(jsonlite)
library(tidyverse)

read_fasta <- function(path) {
    fasta <- as.character(ape::read.FASTA(path))
    fasta <- lapply(fasta, toupper)
    fasta
}

process_file <- function(path) {
    fasta <- read_fasta(path)
    aa <- universal_genetic_code()
    gap_dat <- gap_stats(fasta[[1]], fasta[[2]], aa)

    tab <- gap_dat |> count(phase, type, length)

    #complete(phase = 0:2, type = 1:2, fill = list(n = 0L))

    tab <- tab |> add_column(file = fs::path_file(path), .before = 1L)

    tab
}

gap_stats <- function(ancestor, descendant, genetic_code) {
    aa <- genetic_code

    # Construct codon and amino-acid vectors
    cod_a <- array_to_codons_no_gaps(ancestor)
    cod_d <- array_to_codons_no_gaps(descendant)
    aa_a <- aa[cod_a]
    aa_d <- aa[cod_d]

    # Decode Path information
    # 1 = match; 2 = deletion; 3 = insertion
    aln <- aln_path(ancestor, descendant)
    state <- 1L + ((aln - 1L) %/% 3L)
    state <- c("M", "D", "I")[state]
    phase <- 1L + (aln - 1L) %% 3L

    # Map alignment column to codon position
    pos_a <- (cumsum(state != "I") - 1) %/% 3 + 1
    pos_d <- (cumsum(state != "D") - 1) %/% 3 + 1
    pos_a <- na_if(pos_a, 0)
    pos_d <- na_if(pos_d, 0)

    # Find the coordinates of gaps
    raln <- rle(state)
    rpos <- c(0L, cumsum(raln$lengths))
    rpos1 <- rpos[-length(rpos)] + 1L
    rpos2 <- rpos[-1]
    gaps <- (raln$values != "M") & (raln$length <= 30)

    gap_dat <- tibble(start = rpos1[gaps], end = rpos2[gaps],
        state = raln$values[gaps])
    gap_dat <- gap_dat |> mutate(
        length = end - start + 1L,
        phase = phase[start] - 1L,
        start_cod_a = cod_a[pos_a[start]],
        end_cod_a = cod_a[pos_a[end]],
        start_cod_d = cod_d[pos_d[start]],
        end_cod_d = cod_d[pos_d[end]]
    )
    # add amino acid translations
    gap_dat <- gap_dat |> mutate(
        across(start_cod_a:end_cod_d, \(x) aa[x],
            .names = "{str_replace(.col, 'cod', 'aa')}"))
    gap_dat <- gap_dat |> mutate(type = 1L + (phase != 0 &
        (start_aa_a  != start_aa_d & end_aa_a != end_aa_d)))

    # gap_dat <- gap_dat |> mutate(
    #     phase = factor(phase, 0:2),
    #     type = factor(type, 1:2)
    # )

    gap_dat
}


# encodes an alignment including both phase and alignment state
aln_path <- function(a, d) {
    path <- rep(1:3, length.out = length(a))
    path[d == "-"] <- path[d == "-"] + 3L # deletions
    path[a == "-"] <- path[a == "-"] + 6L # insertions
    path
}

aln_stats <- function(aln) {
    # trim out gaps longer than 30 nucleotides
    raln <- rle( 1 + ((aln-1) %/% 3))
    mask <- (raln$values != 1) & (raln$length > 30)
    raln$values[mask] <- NA
    mask <- inverse.rle(raln)
    aln[is.na(mask)] <- NA

    from <- factor(c(3, aln), levels=1:9)
    to <- factor(c(aln, 1), levels=1:9)
    table(from, to)
}

# g0 <- sum(dat[rbind(c(3 ,4), c(3, 7), c(9, 4))])
# g1 <- sum(dat[rbind(c(1, 5), c(1, 8), c(7, 5))])
# g2 <- sum(dat[rbind(c(2, 6), c(2, 9), c(8, 6))])

array_to_codons_no_gaps <- function(x) {
    x <- x[x != "-"]
    s <- seq(1, length(x), 3)
    str_c(x[s], x[s+1], x[s+2])
}

universal_genetic_code <- function(remove_stops = FALSE) {
    # genetic code in TCGA order
    aa    <- "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    base1 <- "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
    base2 <- "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
    base3 <- "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

    aa <- split_sequence_1(aa)
    base1 <- split_sequence_1(base1)
    base2 <- split_sequence_1(base2)
    base3 <- split_sequence_1(base3)

    names(aa) <- str_c(base1, base2, base3)
    
    if(isTRUE(remove_stops)) {
        aa <- aa[aa != "*"]
    }

    # return code in ACGT order
    aa[order(names(aa))]
}

split_sequence <- function(x) {
    strsplit(x, character(0L))
}

split_sequence_1 <- function(x) {
    split_sequence(x)[[1]]
}

## MAIN ########################################################################

summarize_aln_main <- function(dir) {
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