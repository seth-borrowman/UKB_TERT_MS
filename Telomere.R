library(tidyverse)
library(magrittr)
library(feather)
library(glmnet)
library(survival)
library(ggpubr)
library(rstatix)
set.seed(369)
ts <- read_csv("TSLength.csv", ,
    col_names = TRUE, na = c("", "NA"), trim_ws = TRUE)
filtered_raw <- read_feather("filtered_raw")
covariates <- read_feather("covariates")
covariates <- covariates %>%
    select(all_of(c("Participant ID", "Body mass index (BMI)",
        "IMD", "Fruit", "Veg", "Sex_Male", "Smoking",
        "Summed MET minutes per week for all activity")))

colnames(ts) %<>% str_remove_all(c("\\|")) %>%
    str_remove("Instance 0") %>%
    trimws()

ts <- na.omit(ts)

filtered_raw <- filtered_raw %>% filter(IID %in% ts$`Participant ID`)
ts <- ts %>% filter(`Participant ID` %in% covariates$`Participant ID`)
covariates <- covariates %>% filter(`Participant ID` %in% ts$`Participant ID`)
covariates <- covariates %>% filter(`Participant ID` %in% filtered_raw$IID)
ts <- ts %>% filter(`Participant ID` %in% covariates$`Participant ID`)

x <- merge(ts[, -2], filtered_raw,
    by.x = "Participant ID", by.y = "IID")
x <- merge(x, covariates, by = "Participant ID")
x <- x[order(x$`Participant ID`), ]
x <- as.matrix(x[, -1])
y <- ts[, -3]
y <- y[order(y$`Participant ID`), ]
y <- as.matrix(y[, -1])
y <- log(y)
shapiro.test(y)

cvmod <- cv.glmnet(x, y, family = "gaussian")

coefs_tel <- matrix(rep(0, ncol(x) + 1))
i <- 1
while (i < 1001) {
    samples <- sample(nrow(x), nrow(x), replace = TRUE)
    bootx <- x[samples, ]
    booty <- y[samples, ]
    mod2 <- glmnet(bootx, booty, family = "gaussian", lambda = cvmod$lambda.min)
    coefs_tel <- cbind(coefs_tel, as.matrix(coef(mod2)))
    i <- i + 1
    print(i)
}

coefs_tel <- coefs_tel[, -1]
coefs_tel <- t(coefs_tel)

# Get means of coefficient predictions
means <- rep(0, ncol(coefs_tel))
for (i in seq_len(ncol(coefs_tel))) {
    means[i] <- mean(coefs_tel[, i])
}

# Get SE predictions from bootstrap data
se <- rep(0, ncol(coefs_tel))
for (i in seq_len(ncol(coefs_tel))) {
    thetadot <- sum(coefs_tel[, i] / nrow(coefs_tel))
    se[i] <- sqrt(sum(((coefs_tel[, i] - thetadot)^2) / (nrow(coefs_tel) - 1)))
}

# Create df for plotting coefficient +/- SE predictions
df_tel <- data.frame(Mean = means,
                 StdEr = se,
                 Predictor = colnames(coefs_tel))

df_tel$sig <- as.factor(ifelse((abs(df_tel$Mean) - (1.96 * df_tel$StdEr)) > 0,
     1, 0))

write_csv(df_tel, "df_tel.csv", append = FALSE)

# Plot coefficient estimates and SE
png("plot_tel.png", height = 1080, width = 1920)
ggplot(df_tel, aes(x = Predictor, y = Mean, color = sig)) +
    geom_point(aes(size = 1.5, shape = sig)) +
    geom_errorbar(aes(ymin = Mean - (1.96 * StdEr),
                  ymax = Mean + (1.96 * StdEr)), width = .4,
                  position = position_dodge(0.05)) +
    theme_classic() +
    theme(text = element_text(size = 30), axis.text.x = element_blank(),
     axis.ticks.x = element_blank(), axis.text.y = element_text(size = 30),
     legend.position = "none") +
    scale_color_viridis_d(begin = 0.2, end = 0.7)
dev.off()


merged <- as.data.frame(cbind(x, y))
merged$`5:1268585:C:T_C_2` <- factor(merged$`5:1268585:C:T_C_2`,
    levels = c(0, 1), labels = c("C/T", "C/C"))

c <- colnames(merged)
c[which(colnames(merged) == "5:1268585:C:T_C_2")] <- "rs140124989"
c[which(colnames(merged) == "5:1280296:T:C_T_2")] <- "rs33959226"
colnames(merged) <- c
merged$`rs33959226` <- factor(merged$`rs33959226`,
    levels = c(0, 1), labels = c("T/C", "T/T"))

ttest989 <- merged %>%
    t_test(`Adjusted T/S ratio` ~ rs140124989) %>%
    add_significance()
ttest989 <- ttest989 %>%
    add_xy_position(x = "rs140124989")

png("plot_tel_rs989.png", height = 1080, width = 1080)
merged %>%
    group_by(`rs140124989`) %>%
    mutate(count = n()) %>%
    mutate(mean = mean(`Adjusted T/S ratio`)) %>%
    ggplot(aes(x = `rs140124989`, y = `Adjusted T/S ratio`)) +
    geom_boxplot() +
    ylab("Log Adjusted T/S Ratio") +
    theme_classic() +
    theme(text = element_text(size = 30), axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)) +
    stat_pvalue_manual(ttest989, label = "t-test, p = {p}",
        vjust = -1, label.size = 8)
dev.off()

ttest226 <- merged %>%
    t_test(`Adjusted T/S ratio` ~ rs33959226) %>%
    add_significance()
ttest226 <- ttest226 %>%
    add_xy_position(x = "rs33959226")

png("plot_tel_rs226.png", height = 1080, width = 1080)
merged %>%
    group_by(`rs33959226`) %>%
    mutate(count = n()) %>%
    mutate(mean = mean(`Adjusted T/S ratio`)) %>%
    ggplot(aes(x = `rs33959226`, y = `Adjusted T/S ratio`)) +
    geom_boxplot(lwd = 2, outlier.size = 3) +
    ylab("Log Adjusted T/S Ratio") +
    theme_classic() +
    theme(text = element_text(size = 50), axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40)) +
    stat_pvalue_manual(ttest226, label = "t-test, p = {p}",
        label.size = 10, bracket.size = 2)
dev.off()

a <- select(merged, all_of(c("rs140124989", "rs33959226")))

save.image(file = "TelomereWorkspace.RData")
