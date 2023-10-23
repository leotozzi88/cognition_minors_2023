setwd('/Users/ltozzi/Dropbox (PanLab)/cognition in minors')

data=read.csv('data/Cognition and DASS and BRISC_with z scores merged_minors1 - Cognition and DASS and BRISC_with z scores merged_minors1.csv')

# Drop PTSD
data=data[!data$Group %in% c('ptsd', 'ptsd_hc'), ]
data$Site=factor(data$Site)

#### Create variables of interest
data$interf_t_namecol=data$zvcrtne-data$zvcrtne2 # Create interference time variable

data$maze_perf=rowMeans(data[, c('zemzcomp_r', 'zemzover_r')], na.rm = T)
data$gng_perf=rowMeans(data[, c('zgngavrt_r', 'zgngfp_r', 'zgngfn_r')], na.rm = T)
data$stroop_perf=rowMeans(data[, c('zvi_sco1', 'zvi_sco2')], na.rm = T)
data$vermem_perf=rowMeans(data[, c('zmemtot14', 'zmemrec6', 'zmemrec7')], na.rm = T)
data$contperf_perf=rowMeans(data[, c('zwmrt_r', 'zwmfp_r','zwmfn_r')], na.rm = T)
data$choicert_perf=rowMeans(data[, c('zch_avrt_r', 'ztapdomn', 'ztapdomsd_r')]) # Don't have variability
data$tapping_perf=rowMeans(data[, c('ztapdomn', 'ztapdomsd_r')], na.rm = T) 
data$attswitch_perf=rowMeans(data[, c('zswoadur1_r', 'zesoadur2_r')], na.rm = T) 
data$digitsp_perf=rowMeans(data[, c('zdigitot', 'zrdigitot')], na.rm = T) 

cog_feats=c('maze_perf', 'gng_perf','stroop_perf', 
            'vermem_perf', 'contperf_perf', 'choicert_perf', 
            'tapping_perf', 'attswitch_perf', 'digitsp_perf')

# Remove outliers
data_noout=data
for (feat in cog_feats){
  # Print N removed and %
  print(length(data_noout[!is.na(abs(data[, feat])) & abs(data[, feat])>5, feat]))
  print(length(data_noout[!is.na(abs(data[, feat])) & abs(data[, feat])>5, feat])/nrow(data_noout)*100)
  
  # Remove value
  data_noout[!is.na(abs(data[, feat])) & abs(data[, feat])>5, feat]=NA
}

# Load the reshape2 package
library(reshape2)

# Melt the data into long format
data_noout$Diagnosis=data_noout$Diagnosis_arm
data_noout[data_noout$Clinical==0, 'Diagnosis']='Controls'
data_noout$Diagnosis=factor(data_noout$Diagnosis, levels = c('Controls', 'ADHD', 'Anorexia','FND','First Onset Psychosis'))
mydata_long <- melt(data_noout[, c('Session','Diagnosis', 'Diagnosis_arm', cog_feats, 'NEduS', 'sessage', 'Gender', 'Group2', 'Site')], id.vars = c("Session", 'Diagnosis', 'Diagnosis_arm', 'NEduS', 'sessage', 'Gender', 'Group2', 'Site'), variable.name = "wn", value.name = "Value")

# Run model
library(lme4)
library(lmerTest)

mydata_long$Diagnosis_arm=factor(mydata_long$Diagnosis_arm)
repeated_measures_model=lmer(data=mydata_long, Value~Diagnosis*wn + Site + (1 | Session) )
anova(repeated_measures_model)

# Compute the EMMs for the interaction within cognitive domain
library(emmeans)
emmeans_interaction <- emmeans(repeated_measures_model, ~ Diagnosis | wn, pbkrtest.limit = 7309) # for exact SE use pbkrtest.limit = 7309
posthoc_tests_diags <- pairs(emmeans_interaction, adjust = "tukey")
emmeans_df_withindomain <- as.data.frame(emmeans_interaction)
emmeans_df_withindomain_es=as.data.frame(eff_size(emmeans_interaction, sigma = sigma(repeated_measures_model), edf = df.residual(repeated_measures_model)))

# Compute the EMMs for the interaction within diagnosis
emmeans_interaction <- emmeans(repeated_measures_model, ~ wn | Diagnosis, pbkrtest.limit = 7309) # for exact SE use pbkrtest.limit = 7309
posthoc_tests_cog <- pairs(emmeans_interaction, adjust = "tukey")
emmeans_df_withindiag <- as.data.frame(emmeans_interaction)
emmeans_df_withindiag_es=as.data.frame(eff_size(emmeans_interaction, sigma = sigma(repeated_measures_model), edf = df.residual(repeated_measures_model)))

#### PLOT OF DIFFERENCES FROM CONTROLS

library(ggplot2)

# Extract mean and SE values from emmeans_interaction
repeated_measures_model=lmer(data=mydata_long, Value~Diagnosis*wn + Site + (1 | Session) )
emmeans_interaction <- emmeans(repeated_measures_model, ~ Diagnosis | wn, pbkrtest.limit = 7309) # for exact SE use pbkrtest.limit = 7309
means_emmeans <- summary(emmeans_interaction, infer = TRUE)

# Create a data frame from means_emmeans
data_means_emmeans <- data.frame(
  Diagnosis = means_emmeans$Diagnosis,
  wn = means_emmeans$wn,
  emmean = means_emmeans$emmean,
  SE = means_emmeans$SE
)

data_means_emmeans$wn=factor(data_means_emmeans$wn, levels=levels(data_means_emmeans$wn), labels = c('Executive function', 'Response Inhibition', 'Cognitive Flexibility', 'Verbal Memory', 'Sustained Attention', 'Decision Speed', 'Psychomotor Function', 'Information Processing Speed', 'Working Memory'))

# Plot grouped bar chart with separate panels
p <- ggplot(data = data_means_emmeans, aes(x = Diagnosis, y = emmean, fill = Diagnosis)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", fill = "lightgrey") +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2, position = position_dodge(width = 0.9)) +
  facet_wrap(~ wn, ncol = 3) +
  labs(x = "", y = "Mean Performance (z-score)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.text = element_text(size = 20, face = "bold"),  # Increase size and bold facet labels
    axis.text = element_text(size = 18),# Increase font size for axis labels
    axis.title.y = element_text(size = 20)
  )

# Set the desired font size
png("plots/grouped_bar_chart.png", width = 2000, height = 2500, res = 150)
print(p)
dev.off()


#### PLOT OF PROFILES OF COGNITION

# Plot grouped bar chart with separate panels
p <- ggplot(data = data_means_emmeans, aes(x = wn, y = emmean, fill = Diagnosis)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", fill = "lightgrey") +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2, position = position_dodge(width = 0.9)) +
  facet_wrap(~ Diagnosis, ncol = 5) +
  labs(x = "", y = "Mean Performance (z-score)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.text = element_text(size = 20, face = "bold"),  # Increase size and bold facet labels
    axis.text = element_text(size = 16),# Increase font size for axis labels
    axis.title.y = element_text(size = 16)
  )

# Set the desired font size
png("plots/grouped_bar_chart_cogdomains.png", width = 2500, height = 1000, res = 150)
print(p)
dev.off()



