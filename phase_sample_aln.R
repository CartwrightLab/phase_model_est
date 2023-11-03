OUTPUT_FOLDER <- "results/sampled_aln/"

# Testing
params_file <- "results/params/ATGC001.final.json"
samples_file <- "results/params/ATGC001.samples.rds"

sample_main <- function(samples_file, params_file) {

    dataset <- fs::path_file(samples_file) |> str_extract("^[^.]+(?=[.])")
    withr::local_seed(digest::digest2int(str_glue("{dataset}-20231101")),
        .rng_kind = "Mersenne-Twister")

    # OLD, Bugged solution
    #withr::local_seed(20231101 + digest::digest2int(dataset),
    #    .rng_kind = "Mersenne-Twister")



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
    mat94 <- mg94(params, uni_code61)

    cli_progress_step("Sampling sequences")
    zz_results <- mutate(alignments,
            zz_score_m = zz_score_codons(codons, params, uni_code61),
            zz_score_g = zz_score_gaps(aln, params) ) |>
        mutate( weight = norm_weight(log(n) + 
            zz_score_m + zz_score_g - score),
        .by = file
    )

    # Flush memory
    alignments <- NULL
    gc()

    # Sample Sequences
    zz_results <- slice_sample(zz_results, n = 1, weight_by = weight, by = file)
    zz_results <- zz_results |> select(file, aln, codons) |>
        mutate(file = fs::path_file(file))

    as_seq <- function(file, aln, codons) {
        a <- aln_to_arrays(codons, aln, mat94)

        list(file = file, ancestor = str_flatten(a$ancestor),
            descendant = str_flatten(a$descendant))
    }
    seqs <- pmap(zz_results, as_seq)

    # Write aligned sequences to output files
    cli_progress_step("Writing output fastas")
    out_dir <- fs::dir_create(fs::path(OUTPUT_FOLDER, dataset))

    walk(seqs, function(x) {
        fasta <- str_glue_data(x, ">a\n{ancestor}\n>d\n{descendant}\n")
        write_file(fasta, fs::path(out_dir, x$file))
    })

    invisible()
}

## MODEL #######################################################################

# Use the ZZ model to score alignments
zz_score_codons <- function(codons, params, aa) {
    mat <- log(mg94(params, aa))
    Pi <- log(mg94_stationary_dist(params, aa))
    zz_sub <- t(t(mat) - Pi)

    ret <- codons |> map_dbl(\(x) sum(zz_sub[x]))
    ret
}

zz_score_gaps <- function(aln, params) {
    mat <- zz_gap_matrix(params)

    ret <- map_dbl(aln, function(x) {
        from <- c(3, x)
        to <- c(x, 1)
        sum(mat[cbind(from, to)])
    })
    ret
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

## UTILITIES ###################################################################

aln_path <- function(a, d) {
    path <- rep(1:3, length.out = length(a))
    path[d == "-"] <- path[d == "-"] + 3 # deletions
    path[a == "-"] <- path[a == "-"] + 6 # insertions
    path
}

# reverse the summary stat function from the IS script
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

norm_weight <- function(x) {
    y <- exp(x - max(x))
    y/sum(y)
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
    "rlang", "expm", "digest"))

options(tidyverse.quiet = TRUE)

library(tidyverse)
library(cli)

if(!is_interactive()) {
    args <- commandArgs(trailingOnly = TRUE)
    sample_main(args[[1]], args[[2]])
}
