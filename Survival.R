library(tidyverse)
library(magrittr)
library(glmnet)
library(survival)
library(feather)
set.seed(369)
covariates <- read_feather("covariates")
filtered_raw <- read_feather("filtered_raw")
y <- read_csv("y.csv", col_names = TRUE)
y2 <- read_csv("y2.csv", col_names = TRUE)
cancerdeaths <- read_feather("cancerdeaths")

# Filter covariates and outcomes by those with genotype info
covariates %<>% filter(`Participant ID` %in% filtered_raw$IID)
y %<>% filter(id %in% filtered_raw$IID)
y2 %<>% filter(id %in% filtered_raw$IID)

# Order all data by participant
covariates <- covariates[order(covariates$`Participant ID`), ]
filtered_raw <- filtered_raw[order(filtered_raw$IID), ]
y <- y[order(y$id), ]
y2 <- y2[order(y2$id), ]

# Combine covariates, genotypes
colnames(covariates)[1] <- "ID"
colnames(filtered_raw)[1] <- "ID"
x <- merge(covariates, filtered_raw, by = "ID")

y_crc <- y
y_crc$status <- y_crc$status * 0
y_crc$status[which(y_crc$id %in% cancerdeaths$`Participant ID`)] <- 1
y_crc <- Surv(y_crc$time, y_crc$status)
y_crc_strat <- y2
y_crc_strat <- stratifySurv(y_crc,
    strata = y_crc_strat$strata)

# Change y to a survival object
y <- Surv(y$time, y$status)

y_asc <- y2[which(y2$strata == 1), ]
x_asc <- x %>%
    filter(ID %in% y_asc$id)
y_tra <- y2[which(y2$strata == 2), ]
x_tra <- x %>%
    filter(ID %in% y_tra$id)
y_des <- y2[which(y2$strata == 3), ]
x_des <- x %>%
    filter(ID %in% y_des$id)
y_oth <- y2[which(y2$strata == 0), ]
x_oth <- x %>%
    filter(ID %in% y_oth$id)
y_asc <- Surv(y_asc$time, y_asc$status)
y_tra <- Surv(y_tra$time, y_tra$status)
y_des <- Surv(y_des$time, y_des$status)
y_oth <- Surv(y_oth$time, y_oth$status)

y2 <- stratifySurv(y, strata = y2$strata)

x <- x %>%
    mutate_if(is.character, as.numeric) %>%
    select(!ID) %>%
    as.matrix()
x_asc <- x_asc[-c(4:14)] %>%
    mutate_if(is.character, as.numeric) %>%
    select(!ID) %>%
    as.matrix()
x_des <- x_des[-c(4:14)] %>%
    mutate_if(is.character, as.numeric) %>%
    select(!ID) %>%
    as.matrix()
x_tra <- x_tra[-c(4:14)] %>%
    mutate_if(is.character, as.numeric) %>%
    select(!ID) %>%
    as.matrix()
x_oth <- x_oth[-c(4:14)] %>%
    mutate_if(is.character, as.numeric) %>%
    select(!ID) %>%
    as.matrix()



save.image(file = "SurvivalWorkspace.RData")


cvmod_new <- cv.glmnet(x, y, family = "cox", trace.it = TRUE)

# Non-parametric bootstrap, storing estimates of coefficients
coefs_new <- matrix(rep(0, ncol(x)))
i <- 1
while (i < 1001) { # Update number of iterations
    samples <- sample(nrow(x), nrow(x), replace = TRUE)
    bootx <- x[samples, ]
    booty <- y[samples, ]
    mod2 <- glmnet(bootx, booty, family = "cox", lambda = cvmod_new$lambda.min)
    coefs_new <- cbind(coefs_new, as.matrix(coef(mod2)))
    i <- i + 1
    print(i)
}

# Clean, transpose coefficient estimates
coefs_new <- coefs_new[, -1]
coefs_new <- t(coefs_new)

# Get means of coefficient predictions
means <- rep(0, ncol(coefs_new))
for (i in seq_len(ncol(coefs_new))) {
    means[i] <- mean(coefs_new[, i])
}

# Get SE predictions from bootstrap data
se <- rep(0, ncol(coefs_new))
for (i in seq_len(ncol(coefs_new))) {
    thetadot <- sum(coefs_new[, i] / nrow(coefs_new))
    se[i] <- sqrt(sum(((coefs_new[, i] - thetadot)^2) / (nrow(coefs_new) - 1)))
}

# Create df for plotting coefficient +/- SE predictions
df_new <- data.frame(Mean = means,
                 StdEr = se,
                 Predictor = colnames(coefs_new))
write_csv(df_new, "df_new.csv", append = FALSE)

# Plot coefficient estimates and SE
df_sig <- df_new
df_sig$sig <- as.factor(ifelse((abs(df_sig$Mean) - (1.96 * df_sig$StdEr)) > 0,
     1, 0))

df_sig$Group <- c("Other", "Other", rep("Type of Cancer", 11), "Other",
    rep("Histology of Tumor", 38), rep("Other", 6),
    rep("Treatment", 27),
    rep("SNP", 195))
df_sig <- df_sig %>% mutate(Group = factor(Group, levels = c("SNP",
    "Type of Cancer", "Histology of Tumor",
    "Treatment", "Other")))

png("GroupedCoefficients.png", height = 1080, width = 1920)
ggplot(df_sig, aes(x = Predictor, y = Mean, color = sig)) +
    geom_point(aes(size = 5, shape = sig)) +
    geom_errorbar(aes(ymin = Mean - (1.96 * StdEr),
                  ymax = Mean + (1.96 * StdEr)), width = .6,
                  position = position_dodge(0.05)) +
    theme_classic() +
    theme(text = element_text(size = 50), axis.text.x = element_blank(),
     axis.ticks.x = element_blank(), axis.text.y = element_text(size = 50),
     legend.position = "none") +
    scale_color_viridis_d(begin = 0.2, end = 0.7) +
    facet_grid(~ Group, scales = "free_x")
dev.off()

write.csv(df_sig, "df_sig.csv")

x2 <- x[, -c(3:13)]

cvmod_strat_new <- cv.glmnet(x2, y2, family = "cox", trace.it = TRUE)

# Non-parametric bootstrap, storing estimates of coefficients
coefs_strat_new <- matrix(rep(0, ncol(x2)))
i <- 1
while (i < 1001) { # Update number of iterations
    samples <- sample(nrow(x2), nrow(x2), replace = TRUE)
    bootx <- x2[samples, ]
    booty <- y2[samples, ]
    mod2 <- glmnet(bootx, booty, family = "cox", lambda = cvmod_strat_new$lambda.1se)
    coefs_strat_new <- cbind(coefs_strat_new, as.matrix(coef(mod2)))
    i <- i + 1
    print(i)
}

# Clean, transpose coefficient estimates
coefs_strat_new <- coefs_strat_new[, -1]
coefs_strat_new <- t(coefs_strat_new)

# Get means of coefficient predictions
means <- rep(0, ncol(coefs_strat_new))
for (i in seq_len(ncol(coefs_strat_new))) {
    means[i] <- mean(coefs_strat_new[, i])
}

# Get SE predictions from bootstrap data
se <- rep(0, ncol(coefs_strat_new))
for (i in seq_len(ncol(coefs_strat_new))) {
    thetadot <- sum(coefs_strat_new[, i] / nrow(coefs_strat_new))
    se[i] <- sqrt(sum(((coefs_strat_new[, i] - thetadot)^2) / (nrow(coefs_strat_new) - 1)))
}

# Create df for plotting coefficient +/- SE predictions
df_strat_new <- data.frame(Mean = means,
                 StdEr = se,
                 Predictor = colnames(coefs_strat_new))
write_csv(df_strat_new, "df_strat_new.csv", append = FALSE)

# Plot coefficient estimates and SE
df_sig_strat <- df_strat_new
df_sig_strat$sig <- as.factor(ifelse((abs(df_sig_strat$Mean) - (1.96 * df_sig_strat$StdEr)) > 0,
     1, 0))

df_sig_strat$Group <- c(rep("Other", 3),
    rep("Histology of Tumor", 38), rep("Other", 6),
    rep("Medication or Surgical Treatment", 27),
    rep("SNP", 195))
df_sig_strat <- df_sig_strat %>% mutate(Group = factor(Group, levels = c("SNP",
    "Histology of Tumor",
    "Medication or Surgical Treatment", "Other")))
write.csv(df_sig_strat, "df_sig_strat.csv")

png("GroupedCoefficientsStrat.png", height = 1080, width = 1920)
ggplot(df_sig_strat, aes(x = Predictor, y = Mean, color = sig)) +
    geom_point(aes(size = 1.5, shape = sig)) +
    geom_errorbar(aes(ymin = Mean - (1.96 * StdEr),
                  ymax = Mean + (1.96 * StdEr)), width = .4,
                  position = position_dodge(0.05)) +
    theme_classic() +
    theme(text = element_text(size = 30), axis.text.x = element_blank(),
     axis.ticks.x = element_blank(), axis.text.y = element_text(size = 30),
     legend.position = "none") +
    scale_color_viridis_d(begin = 0.2, end = 0.7) +
    facet_grid(~ Group, scales = "free_x")
dev.off()

save.image(file = "SurvivalWorkspaceNew.RData")

### Ascending

cvmod_asc <- cv.glmnet(x_asc, y_asc, family = "cox", trace.it = TRUE)

# Non-parametric bootstrap, storing estimates of coefficients
coefs_asc <- matrix(rep(0, ncol(x_asc)))
i <- 1
while (i < 1001) { # Update number of iterations
    samples <- sample(nrow(x_asc), nrow(x_asc), replace = TRUE)
    bootx <- x_asc[samples, ]
    booty <- y_asc[samples, ]
    mod2 <- glmnet(bootx, booty, family = "cox", lambda = cvmod_asc$lambda.min)
    coefs_asc <- cbind(coefs_asc, as.matrix(coef(mod2)))
    i <- i + 1
    print(i)
}

# Clean, transpose coefficient estimates
coefs_asc <- coefs_asc[, -1]
coefs_asc <- t(coefs_asc)

# Get means of coefficient predictions
means <- rep(0, ncol(coefs_asc))
for (i in seq_len(ncol(coefs_asc))) {
    means[i] <- mean(coefs_asc[, i])
}

# Get SE predictions from bootstrap data
se <- rep(0, ncol(coefs_asc))
for (i in seq_len(ncol(coefs_asc))) {
    thetadot <- sum(coefs_asc[, i] / nrow(coefs_asc))
    se[i] <- sqrt(sum(((coefs_asc[, i] - thetadot)^2) / (nrow(coefs_asc) - 1)))
}

# Create df for plotting coefficient +/- SE predictions
df_asc <- data.frame(Mean = means,
                 StdEr = se,
                 Predictor = colnames(coefs_asc))
df_asc$Group <- c(rep("Other", 3),
    rep("Histology of Tumor", 38), rep("Other", 6),
    rep("Medication or Surgical Treatment", 27),
    rep("SNP", 195))
df_asc <- df_asc %>% mutate(Group = factor(Group, levels = c("SNP",
    "Histology of Tumor",
    "Medication or Surgical Treatment", "Other")))
df_asc$sig <- as.factor(ifelse((abs(df_asc$Mean) - (1.96 * df_asc$StdEr)) > 0,
     1, 0))
write_csv(df_asc, "df_asc.csv", append = FALSE)

png("GroupedCoefficientsAsc.png", height = 1080, width = 1920)
ggplot(df_asc, aes(x = Predictor, y = Mean, color = sig)) +
    geom_point(aes(size = 1.5, shape = sig)) +
    geom_errorbar(aes(ymin = Mean - (1.96 * StdEr),
                  ymax = Mean + (1.96 * StdEr)), width = .4,
                  position = position_dodge(0.05)) +
    theme_classic() +
    theme(text = element_text(size = 30), axis.text.x = element_blank(),
     axis.ticks.x = element_blank(), axis.text.y = element_text(size = 30),
     legend.position = "none") +
    scale_color_viridis_d(begin = 0.2, end = 0.7) +
    facet_grid(~ Group, scales = "free_x")
dev.off()

### Transverse

cvmod_tra <- cv.glmnet(x_tra, y_tra, family = "cox", trace.it = TRUE)

# Non-parametric bootstrap, storing estimates of coefficients
coefs_tra <- matrix(rep(0, ncol(x_tra)))
i <- 1
while (i < 1001) { # Update number of iterations
    samples <- sample(nrow(x_tra), nrow(x_tra), replace = TRUE)
    bootx <- x_tra[samples, ]
    booty <- y_tra[samples, ]
    mod2 <- glmnet(bootx, booty, family = "cox", lambda = cvmod_tra$lambda.min)
    coefs_tra <- cbind(coefs_tra, as.matrix(coef(mod2)))
    i <- i + 1
    print(i)
}

# Clean, transpose coefficient estimates
coefs_tra <- coefs_tra[, -1]
coefs_tra <- t(coefs_tra)

# Get means of coefficient predictions
means <- rep(0, ncol(coefs_tra))
for (i in seq_len(ncol(coefs_tra))) {
    means[i] <- mean(coefs_tra[, i])
}

# Get SE predictions from bootstrap data
se <- rep(0, ncol(coefs_tra))
for (i in seq_len(ncol(coefs_tra))) {
    thetadot <- sum(coefs_tra[, i] / nrow(coefs_tra))
    se[i] <- sqrt(sum(((coefs_tra[, i] - thetadot)^2) / (nrow(coefs_tra) - 1)))
}

# Create df for plotting coefficient +/- SE predictions
df_tra <- data.frame(Mean = means,
                 StdEr = se,
                 Predictor = colnames(coefs_tra))
df_tra$Group <- c(rep("Other", 3),
    rep("Histology of Tumor", 38), rep("Other", 6),
    rep("Medication or Surgical Treatment", 27),
    rep("SNP", 195))
df_tra <- df_tra %>% mutate(Group = factor(Group, levels = c("SNP",
    "Histology of Tumor",
    "Medication or Surgical Treatment", "Other")))
df_tra$sig <- as.factor(ifelse((abs(df_tra$Mean) - (1.96 * df_tra$StdEr)) > 0,
     1, 0))
write_csv(df_tra, "df_tra.csv", append = FALSE)

png("GroupedCoefficientstra.png", height = 1080, width = 1920)
ggplot(df_tra, aes(x = Predictor, y = Mean, color = sig)) +
    geom_point(aes(size = 1.5, shape = sig)) +
    geom_errorbar(aes(ymin = Mean - (1.96 * StdEr),
                  ymax = Mean + (1.96 * StdEr)), width = .4,
                  position = position_dodge(0.05)) +
    theme_classic() +
    theme(text = element_text(size = 30), axis.text.x = element_blank(),
     axis.ticks.x = element_blank(), axis.text.y = element_text(size = 30),
     legend.position = "none") +
    scale_color_viridis_d(begin = 0.2, end = 0.7) +
    facet_grid(~ Group, scales = "free_x")
dev.off()

### Descending

cvmod_des <- cv.glmnet(x_des, y_des, family = "cox", trace.it = TRUE)

# Non-parametric bootstrap, storing estimates of coefficients
coefs_des <- matrix(rep(0, ncol(x_des)))
i <- 1
while (i < 1001) { # Update number of iterations
    samples <- sample(nrow(x_des), nrow(x_des), replace = TRUE)
    bootx <- x_des[samples, ]
    booty <- y_des[samples, ]
    mod2 <- glmnet(bootx, booty, family = "cox", lambda = cvmod_des$lambda.min)
    coefs_des <- cbind(coefs_des, as.matrix(coef(mod2)))
    i <- i + 1
    print(i)
}

# Clean, transpose coefficient estimates
coefs_des <- coefs_des[, -1]
coefs_des <- t(coefs_des)

# Get means of coefficient predictions
means <- rep(0, ncol(coefs_des))
for (i in seq_len(ncol(coefs_des))) {
    means[i] <- mean(coefs_des[, i])
}

# Get SE predictions from bootstrap data
se <- rep(0, ncol(coefs_des))
for (i in seq_len(ncol(coefs_des))) {
    thetadot <- sum(coefs_des[, i] / nrow(coefs_des))
    se[i] <- sqrt(sum(((coefs_des[, i] - thetadot)^2) / (nrow(coefs_des) - 1)))
}

# Create df for plotting coefficient +/- SE predictions
df_des <- data.frame(Mean = means,
                 StdEr = se,
                 Predictor = colnames(coefs_des))
df_des$Group <- c(rep("Other", 3),
    rep("Histology of Tumor", 38), rep("Other", 6),
    rep("Medication or Surgical Treatment", 27),
    rep("SNP", 195))
df_des <- df_des %>% mutate(Group = factor(Group, levels = c("SNP",
    "Histology of Tumor",
    "Medication or Surgical Treatment", "Other")))
df_des$sig <- as.factor(ifelse((abs(df_des$Mean) - (1.96 * df_des$StdEr)) > 0,
     1, 0))
write_csv(df_des, "df_des.csv", append = FALSE)

png("GroupedCoefficientsdes.png", height = 1080, width = 1920)
ggplot(df_des, aes(x = Predictor, y = Mean, color = sig)) +
    geom_point(aes(size = 1.5, shape = sig)) +
    geom_errorbar(aes(ymin = Mean - (1.96 * StdEr),
                  ymax = Mean + (1.96 * StdEr)), width = .4,
                  position = position_dodge(0.05)) +
    theme_classic() +
    theme(text = element_text(size = 30), axis.text.x = element_blank(),
     axis.ticks.x = element_blank(), axis.text.y = element_text(size = 30),
     legend.position = "none") +
    scale_color_viridis_d(begin = 0.2, end = 0.7) +
    facet_grid(~ Group, scales = "free_x")
dev.off()

### Other

cvmod_oth <- cv.glmnet(x_oth, y_oth, family = "cox", trace.it = TRUE)
cvmod_oth$glmnet.fit
# Non-parametric bootstrap, storing estimates of coefficients
coefs_oth <- matrix(rep(0, ncol(x_oth)))
i <- 1
while (i < 1001) { # Update number of iterations
    samples <- sample(nrow(x_oth), nrow(x_oth), replace = TRUE)
    bootx <- x_oth[samples, ]
    booty <- y_oth[samples, ]
    mod2 <- glmnet(bootx, booty, family = "cox", lambda = cvmod_oth$lambda.min)
    coefs_oth <- cbind(coefs_oth, as.matrix(coef(mod2)))
    i <- i + 1
    print(i)
}

# Clean, transpose coefficient estimates
coefs_oth <- coefs_oth[, -1]
coefs_oth <- t(coefs_oth)

# Get means of coefficient predictions
means <- rep(0, ncol(coefs_oth))
for (i in seq_len(ncol(coefs_oth))) {
    means[i] <- mean(coefs_oth[, i])
}

# Get SE predictions from bootstrap data
se <- rep(0, ncol(coefs_oth))
for (i in seq_len(ncol(coefs_oth))) {
    thetadot <- sum(coefs_oth[, i] / nrow(coefs_oth))
    se[i] <- sqrt(sum(((coefs_oth[, i] - thetadot)^2) / (nrow(coefs_oth) - 1)))
}

# Create df for plotting coefficient +/- SE predictions
df_oth <- data.frame(Mean = means,
                 StdEr = se,
                 Predictor = colnames(coefs_oth))
df_oth$Group <- c(rep("Other", 3),
    rep("Histology of Tumor", 38), rep("Other", 6),
    rep("Medication or Surgical Treatment", 27),
    rep("SNP", 195))
df_oth <- df_oth %>% mutate(Group = factor(Group, levels = c("SNP",
    "Histology of Tumor",
    "Medication or Surgical Treatment", "Other")))
df_oth$sig <- as.factor(ifelse((abs(df_oth$Mean) - (1.96 * df_oth$StdEr)) > 0,
     1, 0))
write_csv(df_oth, "df_oth.csv", append = FALSE)

png("GroupedCoefficientsoth.png", height = 1080, width = 1920)
ggplot(df_oth, aes(x = Predictor, y = Mean, color = sig)) +
    geom_point(aes(size = 1.5, shape = sig)) +
    geom_errorbar(aes(ymin = Mean - (1.96 * StdEr),
                  ymax = Mean + (1.96 * StdEr)), width = .4,
                  position = position_dodge(0.05)) +
    theme_classic() +
    theme(text = element_text(size = 30), axis.text.x = element_blank(),
     axis.ticks.x = element_blank(), axis.text.y = element_text(size = 30),
     legend.position = "none") +
    scale_color_viridis_d(begin = 0.2, end = 0.7) +
    facet_grid(~ Group, scales = "free_x")
dev.off()

# CRC specific death
cvmod_crc <- cv.glmnet(x, y_crc, family = "cox", trace.it = TRUE)

# Non-parametric bootstrap, storing estimates of coefficients
coefs_crc <- matrix(rep(0, ncol(x)))
i <- 1
while (i < 1001) { # Update number of iterations
    samples <- sample(nrow(x), nrow(x), replace = TRUE)
    bootx <- x[samples, ]
    booty <- y_crc[samples, ]
    mod2 <- glmnet(bootx, booty, family = "cox", lambda = cvmod_crc$lambda.1se)
    coefs_crc <- cbind(coefs_crc, as.matrix(coef(mod2)))
    i <- i + 1
    print(i)
}

# Clean, transpose coefficient estimates
coefs_crc <- coefs_crc[, -1]
coefs_crc <- t(coefs_crc)

# Get means of coefficient predictions
means <- rep(0, ncol(coefs_crc))
for (i in seq_len(ncol(coefs_crc))) {
    means[i] <- mean(coefs_crc[, i])
}

# Get SE predictions from bootstrap data
se <- rep(0, ncol(coefs_crc))
for (i in seq_len(ncol(coefs_crc))) {
    thetadot <- sum(coefs_crc[, i] / nrow(coefs_crc))
    se[i] <- sqrt(sum(((coefs_crc[, i] - thetadot)^2) / (nrow(coefs_crc) - 1)))
}

# Create df for plotting coefficient +/- SE predictions
df_crc <- data.frame(Mean = means,
                 StdEr = se,
                 Predictor = colnames(coefs_crc))
df_crc$sig <- as.factor(ifelse((abs(df_crc$Mean) - (1.96 * df_crc$StdEr)) > 0,
     1, 0))
df_crc$Group <- c("Other", "Other", rep("Type of Cancer", 11), "Other",
    rep("Histology of Tumor", 38), rep("Other", 6),
    rep("Medication or Surgical Treatment", 27),
    rep("SNP", 195))
df_crc <- df_crc %>% mutate(Group = factor(Group, levels = c("SNP",
    "Type of Cancer", "Histology of Tumor",
    "Medication or Surgical Treatment", "Other")))
write_csv(df_crc, "df_crc.csv", append = FALSE)

# Plot coefficient estimates and SE
png("plot_crc.png", height = 1080, width = 1920)
ggplot(df_crc, aes(x = Predictor, y = Mean, color = sig)) +
    geom_point(aes(size = 1.5, shape = sig)) +
    geom_errorbar(aes(ymin = Mean - (1.96 * StdEr),
                  ymax = Mean + (1.96 * StdEr)), width = .4,
                  position = position_dodge(0.05)) +
    theme_classic() +
    theme(text = element_text(size = 30), axis.text.x = element_blank(),
     axis.ticks.x = element_blank(), axis.text.y = element_text(size = 30),
     legend.position = "none") +
    scale_color_viridis_d(begin = 0.2, end = 0.7) +
    facet_grid(~ Group, scales = "free_x")
dev.off()

png("SurvCurvCRC.png", height = 1080, width = 1080)
par(mar=c(5,6,4,1)+.1)
plot(survfit(cvmod_crc, x = x, y = y_crc, s = "lambda.1se"),
    xlab = "Time (years)", lwd = 2, bty = "l", cex.lab = 2, cex.axis = 2,
    ylab = "Survival Probability")
dev.off()

# Stratified
cvmod_crc_strat <- cv.glmnet(x2, y_crc_strat,
    family = "cox", trace.it = TRUE)

# Non-parametric bootstrap, storing estimates of coefficients
coefs_crc_strat <- matrix(rep(0, ncol(x2)))
i <- 1
while (i < 1001) { # Update number of iterations
    samples <- sample(nrow(x2), nrow(x2), replace = TRUE)
    bootx <- x2[samples, ]
    booty <- y_crc_strat[samples, ]
    mod2 <- glmnet(bootx, booty, family = "cox",
        lambda = cvmod_crc_strat$lambda.1se)
    coefs_crc_strat <- cbind(coefs_crc_strat, as.matrix(coef(mod2)))
    i <- i + 1
    print(i)
}

# Clean, transpose coefficient estimates
coefs_crc_strat <- coefs_crc_strat[, -1]
coefs_crc_strat <- t(coefs_crc_strat)

# Get means of coefficient predictions
means <- rep(0, ncol(coefs_crc_strat))
for (i in seq_len(ncol(coefs_crc_strat))) {
    means[i] <- mean(coefs_crc_strat[, i])
}

# Get SE predictions from bootstrap data
se <- rep(0, ncol(coefs_crc_strat))
for (i in seq_len(ncol(coefs_crc_strat))) {
    thetadot <- sum(coefs_crc_strat[, i] / nrow(coefs_crc_strat))
    se[i] <- sqrt(sum(((coefs_crc_strat[, i] - thetadot)^2) / (nrow(coefs_crc_strat) - 1)))
}

# Create df for plotting coefficient +/- SE predictions
df_crc_strat <- data.frame(Mean = means,
                 StdEr = se,
                 Predictor = colnames(coefs_crc_strat))
df_crc_strat$sig <- as.factor(ifelse((abs(df_crc_strat$Mean) - (1.96 * df_crc_strat$StdEr)) > 0,
     1, 0))

df_crc_strat$Group <- c(rep("Other", 3),
    rep("Histology of Tumor", 38), rep("Other", 6),
    rep("Medication or Surgical Treatment", 27),
    rep("SNP", 195))
df_crc_strat <- df_crc_strat %>% mutate(Group = factor(Group, levels = c("SNP",
    "Histology of Tumor",
    "Medication or Surgical Treatment", "Other")))
write_csv(df_crc_strat, "df_crc_strat.csv", append = FALSE)

# Plot coefficient estimates and SE
png("plot_crc_strat.png", height = 1080, width = 1920)
ggplot(df_crc_strat, aes(x = Predictor, y = Mean, color = sig)) +
    geom_point(aes(size = 1.5, shape = sig)) +
    geom_errorbar(aes(ymin = Mean - (1.96 * StdEr),
                  ymax = Mean + (1.96 * StdEr)), width = .4,
                  position = position_dodge(0.05)) +
    theme_classic() +
    theme(text = element_text(size = 30), axis.text.x = element_blank(),
     axis.ticks.x = element_blank(), axis.text.y = element_text(size = 30),
     legend.position = "none") +
    scale_color_viridis_d(begin = 0.2, end = 0.7) +
    facet_grid(~ Group, scales = "free_x")
dev.off()

colors <- viridis::viridis(4)
png("SurvCurvStratCRC.png", height = 1080, width = 1080)
par(mar=c(5,6,4,1)+.1)
plot(survfit(cvmod_crc_strat, x = x2, y = y_crc_strat,
    s = "lambda.1se"), col = colors, lty = c(1:4), bty = "l",
    xlab = "Time (years)", lwd = 3, cex.lab = 2, cex.axis = 2,
    ylab = "Survival Probability")
legend("bottomright", legend = c("Ascending", "Transverse",
                                 "Descending", "Other"),
       col = c(colors[2:4], colors[1]), lty = c(2:4, 1),
       bty = "n", cex = 2, lwd = 2)
dev.off()