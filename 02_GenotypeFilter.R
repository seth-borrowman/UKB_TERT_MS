library(tidyverse)
library(magrittr)
covariates <- read_csv("covariates.csv", col_names = TRUE)
raw <- read_delim("rec_snp1.raw", col_names = TRUE)

# Get rid of pesky apostrophes
c <- colnames(covariates)
for (i in 1:length(c)) {
    c[i] <- str_remove_all(c[i], "`")
}
colnames(covariates) <- c

# Remove unnecessary columns
removecols <- c("FID", "PAT", "MAT", "SEX", "PHENOTYPE")
colstokeep <- c(1, seq(2, ncol(raw) - length(removecols), 2))
raw %<>% select(!all_of(removecols)) %>%
    select(all_of(colstokeep))

# Filter genotypes by those in cohort
filtered_raw <- raw %>% filter(IID %in% covariates$`Participant ID`)
covariates <- covariates %>% filter(`Participant ID` %in% filtered_raw$IID)
rm(raw)

# Replace missing values with column average
for (i in 1:ncol(filtered_raw)) {
    filtered_raw[, i][is.na(filtered_raw[, i])] <- unname(colMeans(as.data.frame(filtered_raw[, i]), na.rm = TRUE))
}

# Find columns that have >1 genotype present; create factor
colstofactor <- c()
for (i in 2:ncol(filtered_raw)) {
    if (min(filtered_raw[, i]) != max(filtered_raw[, i])) {
        colstofactor <- append(colstofactor, as.integer(i))
    }
}
filtered_raw[colstofactor] <- lapply(filtered_raw[colstofactor], as.factor)

# Split factored columns for regression
filtered_raw <- model.matrix(~ . - 1, filtered_raw)

# Write to file
write_csv(as.data.frame(filtered_raw), "filtered_raw.csv", append = FALSE)
