library(tidyverse)
library(magrittr)
library(cobalt)
library(feather)
participant_import <- read_csv("Participant_table.csv",
    col_names = TRUE, na = c("", "NA"), trim_ws = TRUE)
drugs <- read_csv("drugs.csv")
fruitveg <- read_csv("FruitVeg.csv")

# Clean up column names
colnames(participant_import) %<>% str_remove_all(c("\\|")) %>%
    str_remove("Instance 0") %>%
    trimws()

# Set date columns to date type
participant_import$`Date of cancer diagnosis` <- as.Date(participant_import$`Date of cancer diagnosis`,
    "%m/%d/%Y")
participant_import$`Date of death` <- as.Date(participant_import$`Date of death`,
    "%m/%d/%Y")
participant_import$`Date of attending assessment centre` <- as.Date(participant_import$`Date of attending assessment centre`,
    "%m/%d/%Y")
participant_import$`Date lost to follow-up` <- as.Date(participant_import$`Date lost to follow-up`,
    "%m/%d/%Y")

# Be sure cohort is prospective and cancer is diagnosed before death
part_prospective <- subset(participant_import,
        `Date of cancer diagnosis` > `Date of attending assessment centre`)
part_diagbeforedeath <- subset(part_prospective,
        `Date of death` != `Date of cancer diagnosis` | is.na(`Date of death`))
part <- subset(part_diagbeforedeath,
        `Date of attending assessment centre` < as.Date("2010-12-31"))
# Remove those with other cancers
for (i in seq_along(part$`Cancer code self-reported`)) {
    part$`Cancer code self-reported`[i] <- str_remove_all(
        part$`Cancer code self-reported`[i], "\\[\"")
    part$`Cancer code self-reported`[i] <- str_remove_all(
        part$`Cancer code self-reported`[i], "\"\\]")
}
part <- subset(part,
        is.na(`Cancer code self-reported`) |
        `Cancer code self-reported` == "skin cancer" |
        `Cancer code self-reported` == "non-melanoma skin cancer" |
        `Cancer code self-reported` == "large bowel cancer/colorectal cancer" |
        `Cancer code self-reported` == "colon cancer/sigmoid cancer\",\"rectal cancer" |
        `Cancer code self-reported` == "colon cancer/sigmoid cancer")

rm(participant_import)
rm(part_prospective)
rm(part_diagbeforedeath)

# Combine IMD
part$IMD <- coalesce(scale(part$`Index of Multiple Deprivation (England)`),
                     scale(part$`Index of Multiple Deprivation (Scotland)`),
                     scale(part$`Index of Multiple Deprivation (Wales)`))
part$IMD[is.na(part$IMD)] <- unname(colMeans(as.data.frame(part$IMD),
    na.rm = TRUE))

# Group cancer sites
for (i in seq_along(part$`Type of cancer: ICD10`)) {
    part$`Type of cancer: ICD10`[i] <- str_remove_all(
        part$`Type of cancer: ICD10`[i], "\\[\"")
    part$`Type of cancer: ICD10`[i] <- str_remove_all(
        part$`Type of cancer: ICD10`[i], "\"\\]")
}
part$site <- rep(0, nrow(part))
asc <- part %>% filter(`Type of cancer: ICD10` %in% c("C18.0 Caecum",
    "C18.1 Appendix", "C18.2 Ascending colon", "C18.3 Hepatic flexure"))
tra <- part %>% filter(`Type of cancer: ICD10` == "C18.4 Transverse colon")
des <- part %>% filter(`Type of cancer: ICD10` %in% c("C18.5 Splenic flexure",
    "C18.6 Descending colon", "C18.7 Sigmoid colon",
    "C19 Malignant neoplasm of rectosigmoid junction",
    "C20 Malignant neoplasm of rectum"))
part$site[which(part$`Participant ID` %in% asc$`Participant ID`)] <- 1
part$site[which(part$`Participant ID` %in% tra$`Participant ID`)] <- 2
part$site[which(part$`Participant ID` %in% des$`Participant ID`)] <- 3

part$Smoking <-
    part$Smoking %>%
    replace(is.na(.), 0)

# Grams of fruits and veggies
fruitveg %<>% filter(`Participant ID` %in% part$`Participant ID`)
part$Fruit <- fruitveg$FreshFruit + fruitveg$DriedFruit
part$Veg <- fruitveg$CookedVeg + fruitveg$RawVeg

# Any comorbidities
part$Comorbid <- ifelse(part$`Long-standing illness disability or infirmity` == "Yes" |
    part$`Other serious medical condition/disability diagnosed by doctor` == "Yes",
    1, 0)

# Operative procedures
for (i in seq_along(part$`Operative procedures - OPCS4`)) {
    part$`Operative procedures - OPCS4`[i] <- str_remove_all(
        part$`Operative procedures - OPCS4`[i], "\\[\"")
    part$`Operative procedures - OPCS4`[i] <- str_remove_all(
        part$`Operative procedures - OPCS4`[i], "\"\\]")
}
ops <- as.data.frame(matrix(rep(0, 72 * nrow(part)),
    nrow = nrow(part), ncol = 72))
ops[, 72] <- part$`Participant ID`
for (i in 1:71) {
    opcs4 <- paste("H", i, sep = "")
    colnames(ops)[i] <- opcs4
    ops[, i] <- ifelse(str_detect(part$`Operative procedures - OPCS4`,
        opcs4), 1, 0)
}
colstoremove <- c()
for (i in seq_len(ncol(ops))){
    if (min(ops[, i], na.rm = TRUE) == max(ops[, i], na.rm = TRUE)) {
            colstoremove <- append(colstoremove, i)
    }
}
ops <- ops[, -colstoremove]
ops <- select(ops, !all_of(c("H1", "H2", "H3", "H35", "H36", "H42",
    "H47", "H48", "H49", "H51", "H52", "H53", "H54", "H55", "H56",
    "H58", "H59", "H60", "H62", "H66", "H68", "H69", "H70")))
for (i in seq_len(ncol(ops))) {
    ops[which(is.na(ops[, i])), i] <- names(which.max(table(ops[, i])))
}
ops <- ops %>%
    apply(2, as.integer) %>%
    as.data.frame()
ops_combined <- as.data.frame(matrix(rep(0, 6 * nrow(part)),
    nrow = nrow(part), ncol = 6))
colnames(ops_combined) <- c("Participant ID", "diag", "endo", "non", "rem", "oth")
ops_combined$`Participant ID` <- ops$V72
ops_combined$diag <- ops$H22 + ops$H28
ops_combined$endo <- ops$H18 + ops$H20 + ops$H21 + ops$H23 + ops$H24 + ops$H25 +
    ops$H26 + ops$H27 + ops$H37
ops_combined$non <- ops$H12 + ops$H34
ops_combined$rem <- ops$H4 + ops$H5 + ops$H6 + ops$H7 + ops$H10 + ops$H11 +
    ops$H14 + ops$H15 + ops$H29 + ops$H32 + ops$H33
ops_combined$oth <- ops$H13 + ops$H16 + ops$H17 + ops$H19 + ops$H30 + ops$H31 +
    ops$H40 + ops$H41 + ops$H44 + ops$H46
ops_combined %<>%
    mutate(across(all_of(names(ops_combined[,-1])), ~ ifelse(.x > 0, 1, 0)))


# Give means or most frequent value to missing observations
numeric <- unname(sapply(part, is.numeric))
character <- unname(sapply(part, is.character))
for (i in seq_len(ncol(part))) {
    if (names(part[, i]) != "Primary cause of death: ICD10" &&
        names(part[, i]) != "Secondary causes of death: ICD10") {
        if (character[i] == TRUE) {
            part[which(is.na(part[, i])), i] <- names(
                which.max(table(part[, i])))
     } else if (numeric[i] == TRUE) {
            part[which(is.na(part[, i])), i] <- unname(
                colMeans(part[, i], na.rm = TRUE))
        }
    }
}

# Make all character columns a factor
part[sapply(part, is.character)] <- lapply(part[sapply(part, is.character)],
                                       as.factor)

# Create new status column for use in surv
status <- rep(0, length(part$`Participant ID`))
for (i in seq_len(length(part$`Participant ID`))) {
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
length2 <- difftime(part$`Date lost to follow-up`,
                    part$`Date of cancer diagnosis`,
                    units = "days")
length2 <- as.numeric(length2) / 365
length3 <- difftime(as.Date("2021-09-30"), part$`Date of cancer diagnosis`,
                    units = "days")
for (i in seq_len(length(length3))) {
    if (status[i] != 0) {
        length3[i] <- NA
    }
}
length3 <- as.numeric(length3) / 365
time <- coalesce(length1, length2, length3)

# New df for surv, includes ID to filter later
id <- part$`Participant ID`
strata <- part$site
y <- data.frame(id, time, status)
y2 <- data.frame(id, time, status, strata)

# Write surv data df to file
write_csv(y, "y.csv", append = FALSE)
write_csv(y2, "y2.csv", append = FALSE)

# Get cancer deaths
deaths <- c("C18", "C19", "C20", "C21.8", "C26.0", "C80",
            "C97", "D37.4", "D37.5")
specdeath <- c()
for (i in seq_along(deaths)) {
    specdeath <- specdeath %>%
        append(which(str_detect(part$`Primary cause of death: ICD10`,
        deaths[i])))
}

cancerdeaths <- part[specdeath, 1]
write_feather(cancerdeaths, "cancerdeaths")

# Column names with covariate data
covariates <- c("Participant ID", "Sex", "Age at cancer diagnosis",
    "Type of cancer: ICD10", "Smoking", "Histology of cancer tumour",
    "Body mass index (BMI)", "IMD", "site",
    "Summed MET minutes per week for all activity",
     "Fruit", "Veg", "Comorbid")

covariates <- part[covariates]

# Split factored columns for regression
covariates <- splitfactor(covariates, drop.first = TRUE)

### Add drug data
drugs <- select(drugs, c("eid", "bnf_code"))
drugs$bnf_code <- as.factor(drugs$bnf_code)
drugs <- splitfactor(drugs, drop.first = FALSE)
across_cols <- names(drugs)
across_cols <- across_cols[-1]
# Multiple entries for each eid, group by eid
drugs %<>%
  group_by(eid) %>%
  summarise(across(all_of(across_cols), sum)) %>%
  mutate(across(all_of(across_cols), ~ ifelse(.x > 0, 1, 0)))
# Only take those in cohort
eid <- covariates$`Participant ID`[!(covariates$`Participant ID` %in% drugs$eid)]
zeros <- as.data.frame(matrix(rep(0, (length(eid) * length(across_cols))),
    nrow = length(eid), ncol = length(across_cols)))
colnames(zeros) <- across_cols
zeros <- data.frame(eid, zeros)
rm(eid)
drugs <- rbind(drugs, zeros)
# Combine with covariates
covariates <- merge(covariates, drugs, by.x = "Participant ID", by.y = "eid")
# Add operations too
covariates <- merge(covariates, ops_combined, by = "Participant ID")

write_feather(as.data.frame(covariates), "covariates")
