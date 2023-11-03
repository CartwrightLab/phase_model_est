# This script uses importance-sampling to estimate parameters for a model of
# coding sequence evolution.

BURNIN_STEPS <- 3L  # execute at least this many steps
MAX_SAMPLE_STEPS <- 10L # maximum number of steps to sample
MAX_STEPS <- 50L    # terminate if this step is reached
SAMPLES_PER_PAIR <- 100L # coati sample size

#### MAIN FUNCTION ###
#
# fasta_directory = directory used to input sequences. Every file ending in
#  .fasta will be included in the analysis

# testing
fasta_directory <- "./data/filtered_cds/02_MusRat"
bin_path <- "/home/reed/Projects/coati/build/release/src/coati-sample"
results_directory <- "./results/params"

phase_im_main <- function(fasta_directory, results_directory, bin_path) {
    # initial parameters
    params0 <- list(
        time = 0.2,
        pi = c(0.25,0.25,0.25,0.25),
        sigma = c(0.25, 0.5, 0.25, 0.25, 0.5, 0.25),
        omega = 0.2,
        gap_open = c(0.001, 0.001, 0.001),
        avg_gap_length = 9
    )

    # setup output file results_dir/dir_name.pid.step.json
    dir_name <- fs::path_file(fs::path_abs(fasta_directory))
    output_stem <- str_c(dir_name, as.integer(now("UTC")), Sys.getpid(), "json", sep="." )
    output_path_stem <- fs::path(fs::path_abs(results_directory), output_stem)

    output_file <- function(step, ext = "json") {
        if(is.numeric(step)) {
            step <- sprintf("%03d", as.integer(step))
        }
        fs::path_ext_set(output_path_stem, str_c(step, ext, sep="."))
    }

    # find the input files and sort based on basename
    fasta_files <- fs::dir_ls(fs::path_abs(fasta_directory), glob="*.fasta")
    #fasta_files <- fasta_files[1000 + (1:10)]

    basenames <- fs::path_file(fasta_files)
    fasta_files <- fasta_files[order(basenames)]

    # construct seeds based on the basenames of input files
    basename_seed <- create_seed(sort(basenames))

    # setup genetic code tables
    uni_code61 <- universal_genetic_code(remove_stops = TRUE)
    codon_list61 <- names(uni_code61)
    amb_code61 <- ambiguous_codons(names(uni_code61), as_index = TRUE)

    step <- 0
    params <- params0
    results <- NULL
    resample <- TRUE
    converged <- FALSE
    document <- list()
    document[[1]] <- list( dataset = dir_name, params = params,
        step = 0, time = now("UTC"), pid = Sys.getpid(),
        resample = resample, converged = converged)

    jsonlite::write_json(document[1], path = output_file(0L),
        digits = NA )

    repeat {
        step <- step + 1
        cli_h1("Step {step}")

        # set R seed, which we will use to generate per-file seeds
        withr::local_seed(808800 + step, .rng_kind = "Mersenne-Twister")

        coati_args <- est_to_coati_args(params, n = SAMPLES_PER_PAIR)
        prev_params <- params

        mat94 <- mg94(params, uni_code61)
        pi94 <- mg94_stationary_dist(params, uni_code61)
        w94 <- mat94*pi94

        if(resample || is.null(results)) {
            progress <- list(
                clear = FALSE,
                type = "custom",
                show_after = 0,
                current = TRUE,
                format = "{cli::pb_spin} Sampling alignments {cli::pb_percent} | ETA: {cli::pb_eta}",
                format_done = "{.alert-success Sampling alignments {.timestamp {cli::pb_elapsed}}}",
                format_failed = "{.alert-danger Sampling alignments {.timestamp {cli::pb_elapsed}}}"
            )

            # Force garbage collection before utilizing a lot of memory
            zz_results <- NULL
            results <- NULL
            gc()
            # Recalculate results
            results <- map(fasta_files, function(filename) {
                my_args <- coati_args
                my_seed <- runif(1L, min = 1, max = 2147483647)
                my_args$seed = as.integer(c(basename_seed, my_seed, 1))
                
                # sample alignments using coati
                out <- call_coati(bin_path, my_args, filename)

                # process alignments and do additional sampling
                aln <- map2(out$ancestor, out$descendant, aln_path)
                codons <- pmap(out, function(ancestor, descendant, score, n) {
                    fill_gaps(ancestor, descendant, score, n, weights = w94,
                        amb_sets = amb_code61)})

                ret <- tibble(aln, codons) |> unnest(2)
                ret
            }, .progress = progress)

            results <- list_rbind(results, names_to = "file")
        }
        
        cli_progress_step("Calculating weights")
        # Force garbage collection before using a lot of memory
        zz_results <- NULL
        gc()
        zz_results <- mutate(results,
                zz_score_m = zz_score_codons(codons, params, uni_code61),
                zz_score_g = zz_score_gaps(aln, params) ) |>
            mutate( weight = norm_weight(log(n) + 
                zz_score_m + zz_score_g - score),
            .by = file
            )

        cli_progress_step("Calculating summary stats")
        codon_counts <- reduce2(zz_results$codons, zz_results$weight,
            function(a, x, w) {
                a + codon_stats(x, wt = w, nbins = length(a))
            }, .init = rep(0, length(mat94)))
        
        aln_counts <- reduce2( zz_results$aln, zz_results$weight,
            function(a, x, w) {
                a + aln_stats(x, wt = w)
            }, .init = rep(0, 81L))

        aln_hist <- reduce2(  zz_results$aln, zz_results$weight,
            function(a, x, w) {
                b <- aln_stats_length_hist(x, wt = w)
                ab <- rep(0, max(length(a), length(b)))
                ab[1L:length(a)] <- a
                ab[1L:length(b)] <- ab[1L:length(b)] + b
                ab
            }, .init = rep(0, 10L))

        cli_progress_step("Estimating parameters")

        # update parameter estimates
        params <- estimate_mg94(codon_counts, params, uni_code61)
        params <- estimate_gaps(aln_counts, params)
        params <- estimate_gap_ext(aln_hist, params)

        # Create and save report
        document[[step + 1]] <- list( dataset = dir_name, params = params,
            step = step, time = now("UTC"), pid = Sys.getpid(),
            resample = resample, converged = converged)

        cli_progress_step("Checking convergence")
        if(step >= BURNIN_STEPS) {
            # After burn-in, test for convergence of substitution model and
            # freeze resampling if it has converged.
            # if sampling is frozen check convergence of all model parameters
            new_params <- list_c(params)
            old_params <- list_c(prev_params)
            delta <- abs(2*(new_params - old_params) / (new_params + old_params))
            if(resample == TRUE) {
                if(max(delta[1:12]) < 0.01 || step >= MAX_SAMPLE_STEPS) {
                    resample <- FALSE
                    # save the frozen sample for later processing
                    saveRDS(results, file =
                        output_file("samples", ext = "rds"))
                }
            } else {
                if(max(delta) < 0.0001) {
                    converged <- TRUE
                }
            }
        }

        document[[step + 1]]$converged <- converged

        cli_progress_step("Saving to file")
        jsonlite::write_json(document[step + 1], path = output_file(step),
            digits = NA )

        cli_progress_done()

        report_params(params, prev_params)

        if( converged ) {
            break
        }
        if(step == MAX_STEPS) {
            converged <- FALSE
            break
        }
    }

    cli_h1("Final Results")

    if(converged) {
        path <- output_file("final")
        jsonlite::write_json(document, path = path, digits = NA)

        cli_alert_success("Estimation converged. Results written to {path}.")
    } else{
        path <- output_file(step)

        cli_alert_danger("Estimation stopped before convergence. Final parameter estimates are in {path}.")
    }

    invisible(document)
}

report_params <- function(new_params, prev_params) {
    o <- c("time", "pi", "sigma", "omega", "gap_open", "avg_gap_length")
    new_params <- list_c(new_params[o])
    prev_params <- list_c(prev_params[o])
    delta <- new_params - prev_params

    cli_div(theme = list("span.myval" = list("color" = "blue"),
        "span.mydelta" = list("color" = "gray")))

    params <- sprintf("{.myval %0.4g} ({.mydelta %+0.2g})",
        new_params, delta) |> map(format_inline)


    cli_h2("Parameter Estimates")
    cli_ul()
    cli_li("Substitution Parameters")
    cli_dl(c("time" = "{params[1]}",
             "pi" = "{params[2:5]}",
             "sigma" = "{params[6:11]}", 
             "omega" = "{params[12]}" ))
    cli_li("Indel Parameters")
    cli_dl(c("gap_open" = "{params[13:15]}",
             "gap_len" = "{params[16]}" ))
}

# Transition probabilities
# 
# m -> m      (1-g)*(1-g)
# m -> d      (1-g)*g
# m -> i      g
# m -> END    (1-g)*(1-g)
# d -> m      1-e
# d -> d      e
# d -> i      0
# d -> END    1-e
# i -> m      (1-e)*(1-g)
# i -> d      (1-e)*g
# i -> i      e
# i -> END    (1-e)*(1-g)

# Transition probabilities for ANC if we strip deletion columns
# m -> m      1-g
# m -> i      g
# m -> END    1-g
# i -> m      1-e
# i -> i      e
# i -> END    1-e

# Transition probabilities for DES if we strip insertion columns
# m -> m      1-g
# m -> d      g
# m -> END    1-g
# d -> m      1-e
# d -> d      e
# d -> END    1-e

# Considering that we want gaps to be in units of three
#
# every N gets a value of (1-g[i])
# every - gets a value of e^(1/3)
# every block of gaps gets a value of g[i]*(1-e)/(1-g[i])/e
#

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

norm_weight <- function(x) {
    y <- exp(x - max(x))
    y/sum(y)
}

# Estimate nucleotide frequencies based on empirical codon frequencies.
# Correct for stop codons not being observed
estimate_pi <- function(seqs) {
    num0 <- map_dbl(c("A", "C", "G", "T"),
        \(x) sum(str_count(seqs, fixed(x))))
    n <- sum(num0)
    Pi <- num0/n
    for(k in 1:20) {
        tag <- Pi[4]*Pi[1]*Pi[3]
        taa <- Pi[4]*Pi[1]*Pi[1]
        tga <- Pi[4]*Pi[3]*Pi[1]
        coding <- 1-(tag+taa+tga)

        # estimate total number of codons after including missing ones
        nn <- n/(3*coding)

        # adjust frequencies
        num <- num0 + nn * c(
            tag+2*taa+tga,
            0,
            tag+tga,
            tag+taa+tga)
        # update Pi
        Pi <- num/sum(num)
    }
    Pi
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
    # trim out gaps longer than 30 nuceotides
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

# h <- aln_hist[seq(3,18,3)]
# m <- weighted.mean(seq_along(h), wt = h)
# p0 <- 1/m

# f <- function(p) {
#     x <- seq_along(h)-1
#     -sum(h*(dgeom(x, prob = p, log = TRUE) - pgeom(max(x), prob = p, log.p = TRUE)))
# }

# estimate mg94+GTR parameters based on codon pairs
estimate_mg94 <- function(counts, params, aa) {
    # helper functions
    params_to_p <- function(params) {
        c(params$pi[1:3]/params$pi[4], params$sigma[1:5]/params$sigma[6],
        params$time, params$omega)
    }

    p_to_params <- function(p) {
        ret <- list(
            pi = c(p[1:3], 1),
            sigma = c(p[4:8], 1),
            time = p[9],
            omega = p[10]
        )
        # fixup pi
        ret$pi <- ret$pi/sum(ret$pi)
        # No need to normalize sigmas every time, but it is here for completeness
        rho <- matrix(0, 4, 4)
        rho[lower.tri(rho)] <- ret$sigma
        rho <- rho + t(rho)
        q_gtr <- t(rho * ret$pi)
        # normalize rates
        ret$sigma <- ret$sigma / sum(q_gtr * ret$pi)
        ret
    }

    # objective function is the -log-likelihood scaled by amount of data
    func <- function(p) {
        par <- p_to_params(p)
        mat <- mg94(par, aa)
        Pi <- mg94_stationary_dist(par, aa)
        w <- -log(mat*Pi)
        val <- sum(w*counts)/sum(counts)
        val
    }
    p0 <- params_to_p(params)
    # save initial value of objective function
    value <- func(p0)
    # use bounded Nelder-Mead optimization
    p_est <- dfoptim::nmkb(p0, func, lower = 0, upper = Inf,
        control=list(tol = 1e-6, trace = FALSE))
    # update parameters on (partial) convergence
    if(p_est$convergence == 0 || p_est$value < value) {
        params <- list_modify(params, !!!p_to_params(p_est$par))
        value <- p_est$value
    }

    params
}

estimate_gaps <- function(aln_counts, params) {
    # let e' = e^(1/3)
    # let (1-e)' = (1-e)/e^(2/3)
    # Gap length probabilities:
    #   3: (1-e)' * e' * e' = (1-e) * e^0
    #   6: (1-e)' * e'^5    = (1-e) * e^(5/3 - 2/3) = (1-e) * e^1

    # geometric: p*(1-p)^(k-1) for k = 1, 2, ...
    #  mean = 1/p, phat = 1/mean(k)

    m2m <- aln_counts[c(3, 10, 20)] # (1-g)(1-g)
    m2d <- aln_counts[c(30, 37, 47)] # (1-g)g
    m2i <- aln_counts[c(57, 64, 74)] # g
    d2m <- aln_counts[c(6, 13, 23)] # (1-e)'
    d2d <- aln_counts[c(33, 40, 50)] # e'
    i2m <- aln_counts[c(9, 16, 26)] # (1-g)*(1-e)'
    i2d <- aln_counts[c(36, 43, 53)] # g*(1-e)'
    i2i <- aln_counts[c(63, 70, 80)] # e'

    gi <- (m2d + m2i + i2d)
    gj <- (2 * m2m + m2d + i2m)
    g <- gi/(gi + gj)

    g <- pmax(g, .Machine$double.eps^0.5)

    ei <- sum(d2d + i2i - 2*(d2m + i2m + i2d))/3
    ej <- sum(d2m + i2m + i2d)

    e <- ei/(ei+ej)

    params <- list_modify(params, gap_open = g, avg_gap_length = 3/(1-e))

    params
}

# estimate gap extension by estimating the average length of gaps of length
# 30 or less
estimate_gap_ext <- function(aln_hist, params) {
    h <- aln_hist[seq(3, 30, 3)]
    m <- weighted.mean( seq_along(h), h)
    a <- if(is.finite(m) && m > 1) 3*m else 3.01
    params <- list_modify(params, avg_gap_length = a)
    params
}

call_coati <- function(bin_path, args, input_file) {
    args <- imap(args, function(x, y) {
        c(str_c("--", y), as.character(x))
    })
    args <- list_c(args)
    args <- c(args, "--", input_file)

    result <- sys::exec_internal(bin_path, args)

    dat <- jsonlite::parse_json(rawToChar(result$stdout))
    # extract data
    score <- map_dbl(dat, "score")
    anc <- map_chr(dat, list("alignment", 1L))
    des <- map_chr(dat, list("alignment", 2L))

    # count 
    tab <- tibble(ancestor = anc, descendant = des, score = score)
    tab <- tab |> count(ancestor, descendant, score)
    tab <- tab |> mutate(ancestor = split_sequence(ancestor),
        descendant = split_sequence(descendant))
    tab
}

# Converts model estimates to closest coati args
est_to_coati_args <- function(params, n = 1L) {
    #--pi, --sigma --omega --time
    #--gap-open, --gap-extend # linear costs
    #--gap-len=3
    #--model=mar-mg
    #--n=1
    #--seed

    # setup constant parameters
    ret <- list("gap-len" = 3, "model" = "mar-mg", "sample-size" = as.integer(n))
    # copy substitution parameters directly
    ret <- c(ret, params[c("time", "pi", "sigma", "omega")])
    # setup gap parameters
    ret[["gap-open"]] <- mean(params$gap_open)
    ret[["gap-extend"]] <- 1-1/params$avg_gap_length

    ret
}

fill_gaps <- function(x, y, score, n, weights, amb_sets) {
    # replace gaps with ambiguous nucleotides
    x[x == "-"] <- "N"
    y[y == "-"] <- "N"
    x <- array_to_codons(x)
    y <- array_to_codons(y)

    # convert codons to indexes
    xn <- match(x, rownames(weights))
    yn <- match(y, colnames(weights))
    is_amb <- is.na(xn) | is.na(yn)

    if(!any(is_amb)) {
        # nothing is ambiguous so we need to do no fill in
        pos <- 1L + (xn - 1L) + (yn - 1L)*nrow(weights)
        tab <- tibble(codons = list(pos), score = score, n = n)
        return(tab)
    }

    # duplicate sequences as necessary
    xn <- matrix(xn, nrow = length(xn), ncol = n)
    yn <- matrix(yn, nrow = length(yn), ncol = n)

    # Randomly fill in codons with ambiguous nucleotides as an additional
    # sampling step.
    for(i in which(is_amb)) {
        xset <- amb_sets[[ x[i] ]]
        yset <- amb_sets[[ y[i] ]]
        w <- weights[xset, yset, drop = FALSE]
        w <- w / sum(w)
        j <- sample(seq_along(w), size = n, prob = w, replace = TRUE) |>
            arrayInd(dim(w))
        xn[i, ] <- xset[ j[,1] ]
        yn[i, ] <- yset[ j[,2] ]
        score <- score + log(w[j])
    }

    pos <- 1L + (xn - 1L) + (yn - 1L)*nrow(weights)
    codons <- apply(pos, 2, c, simplify = FALSE)

    tab <- tibble(codons = codons, score = score) |> count(codons, score)
    tab
}

# encodes an alignment including both phase and alignment state
aln_path <- function(a, d) {
    path <- rep(1:3, length.out = length(a))
    path[d == "-"] <- path[d == "-"] + 3 # deletions
    path[a == "-"] <- path[a == "-"] + 6 # insertions
    path
}

# hash the seed and return it as 4 signed integers
create_seed <- function(strings) {
    # create hash
    h <- rlang::hash(strings)
    # strtoi can't handle signed numbers, use readBin instead
    s <- seq(1, nchar(h), 2)
    n <- strtoi(str_sub(h, s, s+1), base = 16)
    readBin(as.raw(n), integer(), 4, 4, signed=TRUE)
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
    docs <- phase_im_main(args[1], args[2], args[3])
    if(!docs[[length(docs)]]$converged) {
        quit(save = "no", status = 10)
    }
}

# Importance Sampling reference
#  -- https://www.math.arizona.edu/~tgk/mc/book_chap6.pdf