library(tidyverse)
library(magrittr)
ParticipantImport <- read_csv("Participant_table(1).csv", col_names = TRUE, na = c("", "NA"), trim_ws = TRUE)

# Clean up column names
colnames(ParticipantImport) %<>% str_remove_all(c("\\|")) %>%
    str_remove("Instance 0") %>%
    trimws()

# Be sure cohort is prospective and cancer is diagnosed before autopsy
part_prospective <- subset(ParticipantImport,
                           `Date of cancer diagnosis` > `Date of attending assessment centre`)
part_diagbeforedeath <- subset(part_prospective,
                               `Date of death` != `Date of cancer diagnosis` | is.na(`Date of death`))
part <- subset(part_diagbeforedeath,
               `Date of attending assessment centre` < as.Date("2010-12-31"))

rm(ParticipantImport)
rm(part_prospective)
rm(part_diagbeforedeath)

part$IMD <- coalesce(part$`Index of Multiple Deprivation (England)`,
                     part$`Index of Multiple Deprivation (Scotland)`,
                     part$`Index of Multiple Deprivation (Wales)`)
part$IMD[is.na(part$IMD)] <- unname(colMeans(as.data.frame(part$IMD), na.rm = TRUE))

# Make all character columns a factor
part[sapply(part, is.character)] <- lapply(part[sapply(part, is.character)],
                                       as.factor)

# Create new status column for use in surv
status <- rep(0, length(part$`Participant ID`))
for (i in 1:length(part$`Participant ID`)) {
    if (is.na(part[i, "Date of death"])) {
        if (is.na(part[i, "Date lost to follow-up"])) {
            status[i] <- 0
        }
    } else if (is.na(part[i, "Date of death"]) == FALSE) {
        status[i] <- 1
    } else if (is.na(part[i, "Date lost to follow-up"]) == FALSE) {
        status[i] <- 2
    }
}

# Calculate time to event
length1 <- difftime(part$`Date of death`, part$`Date of cancer diagnosis`,
                    units = "days")
length1 <- as.numeric(length1) / 365
length2 <- difftime(part$`Date lost to follow-up`, part$`Date of cancer diagnosis`,
                    units = "days")
length2 <- as.numeric(length2) / 365
length3 <- difftime(Sys.Date(), part$`Date of cancer diagnosis`,
                    units = "days")
for (i in 1:length(length3)) {
    if (status[i] != 0) {
        length3[i] <- NA
    }
}
length3 <- as.numeric(length3) / 365
time <- coalesce(length1, length2, length3)

# New df for surv, includes ID to filter later
ID <- part$`Participant ID`
y <- data.frame(ID, time, status)

# Write surv data df to file
write_csv(y, "y.csv", append = FALSE)

# Column names with covariate data - add as needed -- Include ID to merge with genotype data
covariates <- c("Participant ID", "Sex", "Age at cancer diagnosis", "Type of cancer: ICD10",
                "Pack years adult smoking as proportion of life span exposed to smoking",
                "Alcohol intake frequency.", "Body mass index (BMI)", "IMD")

covariates <- part[covariates]

covariates$`Pack years adult smoking as proportion of life span exposed to smoking` <-
    covariates$`Pack years adult smoking as proportion of life span exposed to smoking` %>%
    replace(is.na(.), 0)

for (i in 1:ncol(covariates)) {
    if (is.factor(covariates[, i]) == TRUE) {
        covariates[, i][is.na(covariates[, i])] <- names(which.max(table(covariates[, i])))
    } else if (is.numeric(covariates[, i]) == TRUE) {
        covariates[, i][is.na(covariates[, i])] <- unname(as.data.frame(colMeans(covariates[, i], na.rm = TRUE)))
    }
}

# Split factored columns for regression
covariates <- model.matrix(~ . - 1, covariates)

write_csv(as.data.frame(covariates), "covariates.csv", append = FALSE)
