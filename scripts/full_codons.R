# USAGE
#This script is designed to traverse a collection of subdirectories  and perform an analysis #on their FASTA files to categorize hypothetical insertions/deletions.  
# This script can be run from R using the `source()` function.
# Or run from bash using `Rscript --vanilla codon-pairs.R
#
# It requires both tidyverse and ape packages to be install.
# They can be installed with the following command:
#  install.packages(c("tidyverse", "ape"))
#





library(tidyverse)
#List subdirectories of current directory
species <- fs::dir_ls(type="dir")

#library(ape)

species_data <- function(j) {
	#files is a list of everything in the current directory that is a file ending in .fasta
	files <- fs::dir_ls(j, type="file",glob="*.fasta")

#Cartwright's code is implemented here:
# For every file, we are going to read the sequences in the file and construct
# a table of codon pairs.
codons <- files |> map(function(f) {
    # read sequences and convert them into character vectors
    dat <- ape::read.FASTA(f)
    dat <- as.character(dat)
    dat <- lapply(dat, tolower)
    # for every sequence, create a sequence of codons, and then create
    # table where each row represents a codon pair. Codons can be used twice
    # once as the before codon and once as the after codon.
    dat <- dat |> map(function(x) {
       # assumes in phase "-"
       dna <- x[x != "-"]
       p <- seq(1, length(dna), 3)
       cod <- str_c(dna[p+0], dna[p+1], dna[p+2])
       tibble(a = head(cod, -1), b = tail(cod, -1))
    })
    dat |> list_rbind()
}, .progress=TRUE)

# Codons is now a data.frame where columns `a` and `b` are adjacent codons
codons <- codons |> list_rbind()

# Construct a named vector that maps codon -> amino acid
aa <- "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
base1  <- "ttttttttttttttttccccccccccccccccaaaaaaaaaaaaaaaagggggggggggggggg"
base2  <- "ttttccccaaaaggggttttccccaaaaggggttttccccaaaaggggttttccccaaaagggg"
base3  <- "tcagtcagtcagtcagtcagtcagtcagtcagtcagtcagtcagtcagtcagtcagtcagtcag"

aa <- aa |> str_split_1("")
base1 <- base1 |> str_split_1("")
base2 <- base2 |> str_split_1("")
base3 <- base3 |> str_split_1("")

names(aa) <- str_glue("{base1}{base2}{base3}")

# Construct a vector that counts every codon pair in the species pair
dat <- codons |> count(a, b)

# Add columns representing what happens if there is a phase 0, phase 1, or
# phase 2 indel
dat <- dat |> mutate(
    f0 = b,
    f1 = str_c(str_sub(a, 1, 1), str_sub(b, 2, 3)),
    f2 = str_c(str_sub(a, 1, 2), str_sub(b, 3, 3))
    )
# Add columns representing the amino acids of the original codons and the
# mutant ones
dat <- dat |> mutate(AA_a = aa[a], AA_b = aa[b],
    AA_f0 = aa[f0], AA_f1 = aa[f1], AA_f2 = aa[f2])

# identify if each indel is of type I or type II.
dat <- dat |> mutate(
    type_f0 = if_else(AA_a == AA_f0 | AA_b == AA_f0, 1, 2),
    type_f1 = if_else(AA_a == AA_f1 | AA_b == AA_f1, 1, 2),
    type_f2 = if_else(AA_a == AA_f2 | AA_b == AA_f2, 1, 2))

# Construct a table that counts how often each phase produces Type I and
# Type II indels.
tab <- dat |>
    filter(!is.na(AA_a) & !is.na(AA_b) & AA_a != "*" & AA_b != "*") |>
    pivot_longer(type_f0:type_f2) |>
    count(name, value, wt = n) |>
    complete(name, value, fill=list(n = 0)) |>
    group_by(name) |> mutate(f = n/sum(n))



}
#end of species_data function

#map species_data to all directories.
	dat <- species  |>
	map(species_data, .progress=TRUE)
	
	#Create dataframe, add column specifying species directory.
	tab <- list_rbind(dat, names_to = "species")
    #write to a file
    #write.csv(tab, "/Users/jacobkebe/Downloads/phase_type_data.csv", row.names = FALSE)

    
 
    
    

