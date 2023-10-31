# This script uses a set of sampled alignments, and parameters to estimate the
# best alignment per sequence pair per species pair

OUTPUT_FOLDER <- "results/best_aln/"

# Testing
params_file <- "results/params/ATGC001.final.json"
samples_file <- "results/params/ATGC001.samples.rds"

best_aln_main <- function(samples_file, params_file) {
    
    cli_progress_step("Reading alignment samples")

    alignments <- readRDS(samples_file)

    cli_progress_step("Reading parameters")

    params <- jsonlite::fromJSON(params_file, simplifyVector = TRUE, 
        simplifyMatrix = FALSE, simplifyDataFrame = FALSE,
        flatten = TRUE)
    params <- params[[length(params)]]
    stopifnot(params$converged)
    params <- params$params

    # setup genetic code tables
    uni_code61 <- universal_genetic_code(remove_stops = TRUE)

    # codon - codon substitution table (includes N)
    w94 <- amb_mg94(params, uni_code61)

    # gap matrix
    wgap <- zz_gap_matrix(params)


    # During the sampling process, gap codons were filled in with nucleotides.
    # We will need to reverse this process.
    cli_progress_step("Generating sequences")
    as_seq <- function(file, aln, codons, score, n) {
        a <- aln_to_arrays(codons, aln, w94)

        # flatten here for better memory consumption (I think. ?)
        tibble(file = file, ancestor = str_flatten(a$ancestor),
            descendant = str_flatten(a$descendant))
    }
    seqs <- pmap(alignments, as_seq)

    # unique alignments
    seqs <- list_rbind(seqs) |> distinct()

    # Score alignments and find the best one per file
    cli_progress_step("Rescoring alignments")
    do_score <- function(anc, des) {
        zz_score(anc, des, w94, wgap)
    }
    seqs <- seqs |> mutate(weight = map2_dbl(ancestor, descendant, do_score),
        file = fs::path_file(file))

    cli_progress_step("Finding best alignments")
    best_aln <- seqs |> slice_max(weight, by = file, with_ties = FALSE)
    best_aln <- transpose(best_aln)


    # Write aligned sequences to output files
    cli_progress_step("Writing output fastas")
    stem <- fs::path_file(samples_file) |> str_extract("^[^.]+(?=[.])")
    out_dir <- fs::dir_create(fs::path(OUTPUT_FOLDER, stem))

    walk(best_aln, function(x) {
        fasta <- str_glue_data(x, ">a\n{ancestor}\n>d\n{descendant}\n")
        write_file(fasta, fs::path(out_dir, x$file))
    })

    invisible()
}

## ZZ MODEL ####################################################################


zz_score <- function(ancestor, descendant, sub_mat, gap_mat) {
    # log-transform substitution matrix
    sub_mat <- log(sub_mat)

    ## score codon patterns

    # sequence to codons
    seq_to_codons <- function(x) {
        x <- split_sequence_1(x)
        x[x == "-"] <- "N"
        array_to_codons(x)
    }

    codons <- rownames(sub_mat)
    stopifnot(codons == colnames(sub_mat))

    anc_s <- seq_to_codons(ancestor)
    des_s <- seq_to_codons(descendant)

    nuc_weight <- sum(sub_mat[cbind(anc_s, des_s)])

    ## score gap patterns
    aln <- aln_path(split_sequence_1(ancestor), split_sequence_1(descendant))
    from <- c(3, aln)
    to <- c(aln, 1)
    gap_weight <- sum(gap_mat[cbind(from, to)])

    ## output
    nuc_weight + gap_weight
}

zz_gap_matrix <- function(params) {
    # phased transition scores
    e <- 3/params$avg_gap_length

    # Adjust gap extension and gap opening scores to score gaps in units of 3
    gap_ext <- log(e) / 3
    gap_end <- log1p(-e)
    gap_open <- log(params$gap_open) - 2*(gap_ext)
    no_gap <- log1p(-params$gap_open)

    mat <- matrix(NA_real_, 9, 9)

    # m -> m     (1-g)*(1-g)
    mat[rbind(c(3, 1), c(1, 2), c(2, 3))] <- 2*no_gap
    # m -> d     (1-g)*g
    mat[rbind(c(3, 4), c(1, 5), c(2, 6))] <- no_gap + gap_open
    # m -> i     g
    mat[rbind(c(3, 7), c(1, 8), c(2, 9))] <- gap_open

    # d -> m     1-e
    mat[rbind(c(6, 1), c(4, 2), c(5, 3))] <- gap_end
    # d -> d     e
    mat[rbind(c(6, 4), c(4, 5), c(5, 6))] <- gap_ext
    # d -> i     0
    mat[rbind(c(6, 7), c(4, 8), c(5, 9))] <- NA_real_

    # i -> m     (1-e)*(1-g)
    mat[rbind(c(9, 1), c(7, 2), c(8, 3))] <- gap_end + no_gap
    # i -> d     (1-e)*g
    mat[rbind(c(9, 4), c(7, 5), c(8, 6))] <- gap_end + gap_open
    # i -> i     e
    mat[rbind(c(9, 7), c(7, 8), c(8, 9))] <- gap_ext

    mat
}

## Utility Functions ###########################################################

# reverse the summary stat function from the sampling script
aln_to_arrays <- function(codons, aln, weights) {
    ind <- arrayInd(codons, dim(weights))

    anc <- rownames(weights)[ind[,1]]
    dec <- colnames(weights)[ind[,2]]
    anc <- split_sequence_1(str_flatten(anc))
    dec <- split_sequence_1(str_flatten(dec))
    aln <- (aln-1) %/% 3

    stopifnot(length(aln) == length(anc))

    dec[aln == 1] <- "-"
    anc[aln == 2] <- "-"

    list(ancestor = anc, descendant = dec)
}

# encodes an alignment including both phase and alignment state
aln_path <- function(a, d) {
    path <- rep(1:3, length.out = length(a))
    path[d == "-"] <- path[d == "-"] + 3 # deletions
    path[a == "-"] <- path[a == "-"] + 6 # insertions
    path
}

## Model Functions #############################################################

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

ambiguous_codons <- function(codons, as_index = FALSE) {
    n1 <- rep(c(FALSE, TRUE), times = 1, each = 4)
    n2 <- rep(c(FALSE, TRUE), times = 2, each = 2)
    n3 <- rep(c(FALSE, TRUE), times = 4, each = 1)

    amb <- split_sequence(set_names(codons)) |> map(function(x) {
        x <- rep(x, each = 8)
        x[c(n1, n2, n3)] <- "N"
        str_c(x[1:8], x[9:16], x[17:24])
    })
    amb <- amb |> map(as_tibble) |>
        list_rbind(names_to="codon")
    amb <- split(amb$codon, amb$value)
    if(isTRUE(as_index)) {
        amb <- map(amb, \(x) match(x, codons))
    }
    amb
}

amb_mg94 <- function(params, genetic_code) {
    codon_list <- names(genetic_code)
    amb_code <- ambiguous_codons(codon_list, as_index = TRUE)
    amb_codons <- names(amb_code)

    # calculate models
    mat94 <- mg94(params, genetic_code)
    pi94 <- mg94_stationary_dist(params, genetic_code)
    w94 <- mat94*pi94

    # construct output matrix
    mat <- matrix(NA_real_, length(amb_codons), length(amb_codons))
    colnames(mat) <- amb_codons
    rownames(mat) <- amb_codons

    for(x in seq_along(amb_codons)) {
        for(y in seq_along(amb_codons)) {
            xx <- amb_code[[ x ]]
            yy <- amb_code[[ y ]]
            mat[x, y] <- sum(w94[xx, yy])
        }
    }
    mat
}

# Qij \propto rho[ic,jc] phi[jc] if A
#             omega * rho[ic,jc] phi[jc] if B
#             0 if otherwise
# where
#   A: i and j are synonymous and differ only at codon position c
#   B: i and j are nonsynonymous and differ only at codon position c

# construct a MG94 matrix based on parameters
mg94 <- function(params, aa) {
    codons <- names(aa)
    tab <- expand_grid(x = codons, y = codons)
    tab <- tab |> mutate(syn = aa[x] == aa[y])

    # construct gtr matrix
    rho <- matrix(0, 4, 4)
    rho[lower.tri(rho)] <- params$sigma
    rho <- rho + t(rho)
    q_gtr <- t(rho * params$pi)
    # normalize rates
    q_gtr <- q_gtr/sum(q_gtr * params$pi)

    # identify which position differs or return NA
    find_pos <- function(x, y) {
        b1 <- str_sub(x, 1, 1) != str_sub(y, 1, 1)
        b2 <- str_sub(x, 2, 2) != str_sub(y, 2, 2)
        b3 <- str_sub(x, 3, 3) != str_sub(y, 3, 3)

        if_else(b1 + b2 + b3 == 1, b1 * 1L + b2 * 2L + b3 * 3L,
            NA_integer_)
    }
    NUC <- c("A" = 1, "C" = 2, "G" = 3, "T" = 4)
    tab <- tab |> mutate(pos = find_pos(x, y))
    tab <- tab |> mutate(x_nuc = NUC[str_sub(x, pos, pos)],
        y_nuc = NUC[str_sub(y, pos, pos)])

    array_ind <- cbind(tab$x_nuc, tab$y_nuc)

    # mg94 Q matrix
    tab <- tab |> mutate(Q = if_else(syn, 1, params$omega)*q_gtr[array_ind])
    tab <- tab |> mutate(Q = coalesce(Q, 0))
    Q <- matrix(tab$Q, length(codons), byrow=TRUE)
    diag(Q) <- -rowSums(Q)

    mat <- expm::expm(Q*params$time)
    rownames(mat) <- codons
    colnames(mat) <- codons
    mat
}



mg94_stationary_dist <- function(params, aa) {
    codons <- names(aa)

    NUC <- c("A" = 1, "C" = 2, "G" = 3, "T" = 4)
    x1 <- NUC[str_sub(codons, 1, 1)]
    x2 <- NUC[str_sub(codons, 2, 2)]
    x3 <- NUC[str_sub(codons, 3, 3)]

    Pi <- params$pi[x1] * params$pi[x2] * params$pi[x3]
    Pi/sum(Pi)
}

codon_stats <- function(codons, wt = 1, nbins = 61*61) {
    tabulate(codons, nbins = nbins)*wt
}

aln_stats <- function(aln, wt = 1) {
    # trim out gaps longer than 12 nuceotides
    raln <- rle( 1 + ((aln-1) %/% 3))
    mask <- (raln$values != 1) & (raln$length > 30)
    raln$values[mask] <- NA
    mask <- rep(raln$values, times = raln$length)
    aln[is.na(mask)] <- NA

    from <- c(3, aln) - 1
    to <- c(aln, 1) - 1
    tabulate(1 + from + 9*to, nbins = 81)*wt
}

aln_stats_length_hist <- function(aln, wt = 1) {
    aln <- c(3, aln, 1)
    aln <- 1 + ((aln-1) %/% 3)
    r <- rle(aln)
    tabulate(r$lengths[r$values > 1])*wt
}

array_to_codons <- function(x) {
    s <- seq(1, length(x), 3)
    str_c(x[s], x[s+1], x[s+2])
}

split_sequence <- function(x) {
    # base R is faster than stringr here
    strsplit(x, character(0L))
}

split_sequence_1 <- function(x) {
    split_sequence(x)[[1]]
}

## HELPERS #####################################################################

# Source: https://github.com/r-lib/cpp11/blob/main/R/utils.R
# Copyright (c) 2020 RStudio (MIT License)

stop_unless_installed <- function(pkgs) {
  has_pkg <- logical(length(pkgs))
  for (i in seq_along(pkgs)) {
    has_pkg[[i]] <- requireNamespace(pkgs[[i]], quietly = TRUE)
  }
  if (any(!has_pkg)) {
    msg <- sprintf(
      "The %s package(s) are required for this functionality",
      paste(pkgs[!has_pkg], collapse = ", ")
    )

    if (is_interactive()) {
      ans <- readline(paste(c(msg, "Would you like to install them? (Y/N) "), collapse = "\n"))
      if (tolower(ans) == "y") {
        utils::install.packages(pkgs[!has_pkg])
        stop_unless_installed(pkgs)
        return()
      }
    }

    stop(msg, call. = FALSE)
  }
}

is_interactive <- function() {
    opt <- getOption("rlang_interactive", NULL)
    if (!is.null(opt)) {
        return(opt)
    }
    if (isTRUE(getOption("knitr.in.progress"))) {
        return(FALSE)
    }
    if (isTRUE(getOption("rstudio.notebook.executing"))) {
        return(FALSE)
    }
    if (identical(Sys.getenv("TESTTHAT"), "true")) {
        return(FALSE)
    }
    interactive()
}

## MAIN ACTION #################################################################

stop_unless_installed(c("tidyverse", "cli", "fs", "withr", "jsonlite", "sys", 
    "rlang", "dfoptim", "expm"))

options(tidyverse.quiet = TRUE)

library(tidyverse)
library(cli)

if(!is_interactive()) {
    args <- commandArgs(trailingOnly = TRUE)
    best_aln_main(args[1], args[2])
}