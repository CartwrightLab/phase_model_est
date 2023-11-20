options(tidyverse.quiet = TRUE, readr.show_col_types = FALSE)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

dat <- read_csv(args[[1]])

tab <- dat |> count(species, length, wt = n)

tab_phase <- dat |> count(species, length, phase, wt = n) |>
    pivot_wider(names_from = "phase", names_prefix = "n_phase_",
        values_from = "n", values_fill = 0)

tab_type <- dat |> count(species, length, type, wt = n) |>
    pivot_wider(names_from = "type", names_prefix = "n_type_",
        values_from = "n", values_fill = 0)

tab <- tab |> inner_join(tab_phase) |>
    inner_join(tab_type)

cat(format_csv(tab))
