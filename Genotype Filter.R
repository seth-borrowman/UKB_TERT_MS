library(tidyverse)
library(magrittr)
library(feather)
library(zoo)
library(cobalt)
covariates <- read_feather("covariates")
raw <- read_delim("snps.raw", col_names = TRUE)

# Get rid of pesky apostrophes
c <- colnames(covariates)
for (i in seq_len(length(c))) {
    c[i] <- str_remove_all(c[i], "`")
}
colnames(covariates) <- c

# Remove unnecessary columns
removecols <- c("FID", "PAT", "MAT", "SEX", "PHENOTYPE")
raw %<>% select(!all_of(removecols))

# Filter genotypes by those in cohort
filtered_raw <- raw %>% filter(IID %in% covariates$`Participant ID`)
covariates <- covariates %>% filter(`Participant ID` %in% filtered_raw$IID)
rm(raw)

# Replace missing values with column average
filtered_raw <- na.aggregate(filtered_raw, FUN = median)

# Find columns that have >1 genotype present; create factor
colstofactor <- c()
for (i in 2:ncol(filtered_raw)) {
    if (min(filtered_raw[, i]) != max(filtered_raw[, i])) {
        colstofactor <- append(colstofactor, as.integer(i))
    }
}
filtered_raw %<>% select(c(`IID`, all_of(colstofactor)))
a <- 2:ncol(filtered_raw)
filtered_raw[, a] <- lapply(filtered_raw[, a], as.factor)

# Split factored columns for regression
filtered_raw <- splitfactor(filtered_raw, drop.first = TRUE)

# Write to file
write_feather(as.data.frame(filtered_raw), "filtered_raw")

5:1268585
filtered_raw <- filtered_raw[-2778, ]
table(filtered_raw$`5:1268585:C:T_C`)
