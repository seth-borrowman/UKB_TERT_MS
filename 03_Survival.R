library(tidyverse)
library(magrittr)
library(glmnet)
library(survival)
covariates <- read_csv("covariates.csv", col_names = TRUE)
filtered_raw <- read_csv("filtered_raw.csv", col_names = TRUE)
y <- read_csv("y.csv", col_names = TRUE)

# Get rid of pesky apostrophes
c <- colnames(covariates)
for (i in 1:length(c)) {
    c[i] <- str_remove_all(c[i], "`")
}
colnames(covariates) <- c
c <- colnames(filtered_raw)
for (i in 1:length(c)) {
    c[i] <- str_remove_all(c[i], "`")
}
colnames(filtered_raw) <- c
c <- colnames(y)
for (i in 1:length(c)) {
    c[i] <- str_remove_all(c[i], "`")
}
colnames(y) <- c


# Filter covariates and outcomes by those with genotype info
covariates %<>% filter(`Participant ID` %in% filtered_raw$IID)
y %<>% filter(ID %in% filtered_raw$IID)

# Order all data by participant
covariates <- covariates[order(covariates$`Participant ID`), ]
filtered_raw <- filtered_raw[order(filtered_raw$IID), ]
y <- y[order(y$ID), ]

# Change y to a survival object
y <- Surv(y$time, y$status)

# Combine covariates, genotypes
colnames(covariates)[1] <- "ID"
colnames(filtered_raw)[1] <- "ID"
x <- merge(covariates, filtered_raw, by = "ID")
x <- as.matrix(x)

# 10-fold CV for double checking model works
# Update nfolds
cvmod <- cv.glmnet(x, y, family = "cox", type.measure = "C", trace.it = TRUE, nfolds = 4)

# Non-parametric bootstrap, storing estimates of coefficients
coefs <- matrix(rep(0, ncol(x)))
i <- 1
while (i < 1000) { # Update number of iterations
    samples <- sample(nrow(x), nrow(x), replace = TRUE)
    bootx <- x[samples, ]
    booty <- y[samples, ]
    #mod1 <- cv.glmnet(bootx, booty, family = "cox", type.measure = "C")
    mod2 <- glmnet(bootx, booty, family = "cox", lambda = cvmod$lambda.1se)
    coefs <- cbind(coefs, as.matrix(coef(mod2)))
    i <- i + 1
}

# Clean, transpose coefficient estimates
coefs <- coefs[, -1]
coefs <- t(coefs)

# Get means of coefficient predictions
means <- rep(0, ncol(coefs))
for (i in 1:ncol(coefs)) {
    means[i] <- mean(coefs[, i])
}

# Get SE predictions from bootstrap data
se <- rep(0, ncol(coefs))
for (i in 1:ncol(coefs)) {
    thetadot <- sum(coefs[, i] / nrow(coefs))
    se[i] <- sqrt(sum(((coefs[, i] - thetadot)^2) / (nrow(coefs) - 1)))
}

# Create df for plotting coefficient +/- SE predictions
df <- data.frame(Mean = means,
                 StdEr = se,
                 Predictor = colnames(coefs))

# Plot coefficient estimates and SE
png("plot.png", height = 1080, width = 1920)
ggplot(df, aes(x = Predictor, y = Mean, label = Predictor)) +
    geom_point() +
    geom_errorbar(aes(ymin = Mean - StdEr, ymax = Mean + StdEr), width = .2,
                  position = position_dodge(0.05)) +
    geom_text(aes(label = ifelse((abs(Mean) - StdEr) > 0,
     as.character(Predictor), "")), hjust = 1, vjust = 0, size = 12) +
    theme_classic() +
    theme(text = element_text(size = 30), axis.text.x = element_text(size = 0),
     axis.text.y = element_text(size = 20)) +
    ggtitle("Mean coefficient estimate with 1000 bootstrap SE")
dev.off()
