library(ggplot2)
library(magrittr)
library(dplyr)
library(latex2exp)
library(survival)
library(glmnet)
load("SurvivalWorkspaceNew.RData")
#read_csv("df_sig_strat.csv")

df_sig <- df_new
df_sig$sig <- as.factor(ifelse((abs(df_sig$Mean) - (1.96 * df_sig$StdEr)) > 0,
     1, 0))

df_sig$Group <- c("Other", "Other", rep("Type of Cancer", 11), "Other",
    rep("Histology of Tumor", 38), rep("Other", 6),
    rep("Medication or Surgical Treatment", 27),
    rep("SNP", 195))
df_sig <- df_sig %>% mutate(Group = factor(Group, levels = c("SNP",
    "Type of Cancer", "Histology of Tumor",
    "Medication or Surgical Treatment", "Other")))

png("GroupedCoefficients.png", height = 1080, width = 1920)
ggplot(df_sig, aes(x = Predictor, y = Mean, color = sig)) +
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

write.csv(df_sig, "df_sig.csv")

df_sig_strat <- df_strat_new
df_sig_strat$sig <- as.factor(ifelse((abs(df_sig_strat$Mean) - (1.96 * df_sig_strat$StdEr)) > 0,
     1, 0))

df_sig_strat$Group <- c(rep("Other", 3),
    rep("Histology of Tumor", 38), rep("Other", 5),
    rep("Medication or Surgical Treatment", 65),
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

png("SurvCurv.png", height = 1080, width = 1080)
par(mar=c(5,6,4,1)+.1)
plot(survfit(cvmod_new, x = x_new, y = y_new, s = "lambda.min"),
    xlab = "Time (years)", lwd = 2, bty = "l", cex.lab = 2, cex.axis = 2,
    ylab = "Survival Probability")
dev.off()

colors <- viridis::viridis(4)
png("SurvCurvStrat.png", height = 1080, width = 1080)
par(mar=c(5,6,4,1)+.1)
plot(survfit(cvmod_strat_new, x = x2, y = y2,
    s = "lambda.min"), col = colors, lty = c(1:4), bty = "l",
    xlab = "Time (years)", lwd = 4, cex.lab = 3, cex.axis = 2,
    ylab = "Survival Probability")
legend("bottomright", legend = c("Ascending", "Transverse",
                                 "Descending", "Other"),
       col = c(colors[2:4], colors[1]), lty = c(2:4, 1),
       bty = "n", cex = 2.5, lwd = 3)
dev.off()


causes <- readr::read_csv("Participant_table(1).csv")
causes <- select(causes, all_of((c("Participant ID",
    "Underlying (primary) cause of death: ICD10 | Instance 0",
    "Contributory (secondary) causes of death: ICD10 | Instance 0 | Array 1"))))
causes <- causes %>% filter(`Participant ID` %in% filtered_raw$ID)
blanks <- which(is.na(causes$`Underlying (primary) cause of death: ICD10 | Instance 0`))
causes <- causes[-blanks, ]
readr::write_csv(causes, "causes.csv")

plot(survfit(cvmod_strat_new, x = x2_new, y = y2_new,
    s = "lambda.min"), fun = "cloglog")
