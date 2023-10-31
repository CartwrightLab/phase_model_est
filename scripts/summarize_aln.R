options(tidyverse.quiet = TRUE)
library(jsonlite)
library(tidyverse)

files <- fs::dir_ls(".", type="file", glob="*.fasta")

f <- files[1]

dat <- ape::read.FASTA(f)
dat <- as.character(dat)
dat <- lapply(dat, toupper)
