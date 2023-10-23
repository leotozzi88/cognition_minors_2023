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

# Get the true covariance matrix
cov_matrix=cov(data[, cog_feats], use = 'pairwise.complete.obs')

#### Simulation
set.seed(123)  # for reproducibility

# Initialize some variables
n_runs <- 1000  # Number of simulations, set to 1 for this example
diagnosis_types <- c("Control", "ADHD", "Anorexia", "First Onset Psychosis", "FND")
n_per_group <- c(483, 343, 40, 25, 56)
tests <- c("Executive Function", "Response Inhibition", "Cognitive Flexibility", "Verbal Memory", "Sustained Attention", "Decision Speed", "Psychomotor Function", "Information Processing Speed", "Working Memory")
dysfs=seq(0, 1, by=0.1) # dysfunctions to test
ndysfs=seq(2, 8, by=2) # number of dysfunctional domains


# Initialize dataframe to store power calculations
power_df <- data.frame()

for (ndysf in ndysfs){
  
  # Generate dysfunction patterns for each group
  mask_mat <- matrix(0, nrow = 5, ncol = 9, dimnames = list(diagnosis_types))
  # For each row starting from the second row, sample a sequence of 4 1s and 5 0s
  for (i in 2:5) {
    mask_mat[i, ] <- sample(c(rep(1, ndysf), rep(0, length(tests)-ndysf)))
  }
  
  for (dysf in dysfs){
    
    allrej=c() #Vector to store rejections of null
    
    for (run in 1:n_runs){
      
      print(paste('ndysf', ndysf,'dysf', dysf, 'run', run))
      
      # Create the simulated data frame
      repeated_diagnosis <- rep(diagnosis_types, times = n_per_group)
      subidx=1:length(repeated_diagnosis)
      sim_data=data.frame(id=subidx, diag=repeated_diagnosis)
      
      for (sub in 1:length(subidx)){
        
        # Generate the mean performance of the subject based on their diagnosis
        subperf_means=mask_mat[sim_data[sub, 'diag'], ]*dysf
        
        # Generate a vector of dysfunctions
        perf=as.numeric(rmvnorm(n=1, mean = subperf_means, sigma = cov_matrix))
        
        # Assign to the data frame
        sim_data[sub, tests]=perf
        
      }
      
      # Convert to long format
      library(tidyverse)
      long_sim_data <- sim_data %>% 
        pivot_longer(
          cols = -c(id, diag), # Columns to reshape
          names_to = "tests",  # New column for test names
          values_to = "performance" # New column for test values
        )
      
      # Run model
      library(lme4)
      library(lmerTest)
      repeated_measures_model=lmer(data=long_sim_data, performance~diag*tests + (1 | id) )
      
      # Extract interaction p-value
      pval=anova(repeated_measures_model)[3, 6]
      allrej=append(allrej, (pval<0.05)*1)
      
    }
    
    # Add information to power data frame
    power_df=rbind(power_df, data.frame(ndysf=ndysf, dysf=dysf, power=mean(allrej)))
    
  }
}

write.csv(power_df, 'tables/power_df.csv', row.names=F)

# Create the line plot
power_df$ndysf=factor(power_df$ndysf)
power_df$dysf=factor(power_df$dysf)
ggplot(power_df, aes(x = dysf, y = power, color=ndysf, group = ndysf)) +
  geom_line(size=1) +
  labs(
    x = "Mean Performance Difference vs. Controls for Dysfunctional Cognitive Domains",
    y = "Power for Diagnosis * Cognitive Domain Interaction",
    color = "Number of Dysfunctional Cognitive Domains:"
  ) +
  theme_minimal(base_size = 16)+
  theme(legend.position = "top")




