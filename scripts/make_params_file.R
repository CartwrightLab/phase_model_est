options(tidyverse.quiet = TRUE)
library(jsonlite)
library(tidyverse)

dat <- read_json("stdin", simplifyVector = TRUE, flatten = TRUE)

dat <- as_tibble(dat) |> select(dataset, starts_with("params"))

dat <- dat |> rename_with(str_remove, pattern = "^params.")

dat <- dat |> mutate(across(c(dataset, time, omega,
    avg_gap_length), list_c))

dat <- dat |> unnest_wider(pi, names_sep=".")
dat <- dat |> unnest_wider(sigma, names_sep=".")
dat <- dat |> unnest_wider(gap_open, names_sep=".")

dat <- dat |> mutate(rho.phase0 = -log1p(-dat$gap_open.1)/time,
    rho.phase1 = -log1p(-dat$gap_open.2)/time,
    rho.phase2 = -log1p(-dat$gap_open.3)/time) |>
    select(-starts_with("gap_open"))

dat <- dat |> rename(branch_length = time,
    pi.a = pi.1, pi.c = pi.2, pi.g = pi.3, pi.t = pi.4,
    sigma.ac = sigma.1, sigma.ag = sigma.2, sigma.at = sigma.3,
    sigma.cg = sigma.4, sigma.ct = sigma.5, sigma.gt = sigma.6)

dat <- dat |> relocate(avg_gap_length, .after = -1)

cat(format_csv(dat))
