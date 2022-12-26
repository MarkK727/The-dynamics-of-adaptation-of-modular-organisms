######################## Supplementary figure - standard deviation of WT-GFP brightness measurements ########################

library(ggplot2)
library(ggpubr)
library(gridExtra)
library(stringr)
library(tidyverse)

# Save the nucleotide-genotype to brightness dataset as df_GFP
df_GFP = read.csv("Desktop/GFP Dataset/df.avGFP.csv")
# Save a separate dataframe with or w/o the wildtype data. E.g. df_GFP$Mut_number > 0, if you want to omit WT.
df_GFP_omit_WT = df_GFP[which(df_GFP$Mut_number >= 0),]

# Calculate the standard errors if possible (i.e. uniqueBarcodes > 1). SE = sig/sqrt(n)
replicate_indices = which(df_GFP_omit_WT$uniqueBarcodes > 1, arr.ind=TRUE)
df_GFP_replicates = df_GFP_omit_WT[replicate_indices,]
SE = as.data.frame(df_GFP_replicates$std/sqrt(df_GFP_replicates$uniqueBarcodes))
names(SE) = "SE"
df_GFP_replicates[,9]= SE
df_GFP_replicates = as.data.frame(df_GFP_replicates)

x1 = as.factor(df_GFP_replicates$uniqueBarcodes)
x2 = as.factor(df_GFP_replicates$medianBrightness)
y = as.factor(df_GFP_replicates$SE)

SI_std1 = ggplot(df_GFP_replicates, aes(x=log10(uniqueBarcodes), y=std)) + geom_point(size=2, shape=1) + geom_smooth(method = lm,se=FALSE, fullrange=TRUE) +
  theme_bw(base_size = 10) + xlab("The Number of Replicate Measurements (log10)") + ylab("Standard Deviation")

# The number of unique barcodes vs. Standard error
SI_std2 = ggplot(df_GFP_replicates, aes(x=log10(uniqueBarcodes), y=SE)) + geom_point(size=2, shape=1) + geom_smooth(method = lm,se=FALSE, fullrange=TRUE) +
  theme_bw(base_size = 10) + xlab("The Number of Replicate Measurements (log10)") + ylab("Standard Error")

SI_std3 = ggplot(df_GFP_replicates, aes(std)) + geom_histogram(aes(y = after_stat(density)), 
                                                          colour = "black", 
                                                          fill = "grey") + 
  theme_bw(base_size = 10) + xlab("Standard deviation") + ylab("density") +
  stat_function(fun = dnorm, args = list(mean = mean(df_GFP_replicates$std), sd = sd(df_GFP_replicates$std)), colour = "blue")
  

######################## Figure 5 - standard deviation of WT-GFP brightness measurements ######################## 


# Benjamini-Hochberg critical value
FDR_critVal = 0.05
bin_parameter = 1

df = df_GFP

# First, for each row containing a mutant strain, add a column that contains the list of neighbor index ( e.g. c(2,45,68) )
# Calculate length of the list
N = length(df[[1]])

Mut_number = df$Mut_number
# Identify all unique n for nth-mutant
n = unique(Mut_number)

# Break the full dataset into sub-dataset of nth mutant classes. e.g. Dataframe "1" contains all 1-mutant from wt. "0" is wt. 
# All names (e.g. 0,1,2,... ) contained in the "Name_List"
t = 1
Name_List = c()

for (n_i in n) {
  nam <- paste("", n_i, sep = "")
  Name_List[t] = nam
  
  Neighbor_Index = which(df$Mut_number == n_i)
  
  dat = df[which(df$Mut_number == n_i),]
  assign(nam, dat)
  
  t = t+1
}

l = sort(as.numeric(Name_List))

p_list = df$pval
p_vec = c()
N_tot = c()
p_adjusted = c()
bFrac = c()
N_benef = c()

k_list = df$k
k_vec = c()
kb_list = c()
kb_mean = c()
kb_std = c()

for ( i in 1:length(p_list) ) {
  # Convert the character string to double (fraction of muational neighbors with higher performing phenotype relative to the focal genotype)
  p_vec = eval(parse(text = p_list[[i]]))
  BH_corr = p.adjust(p_vec, method = "BH", n = length(p_vec))
  p_adjusted[i] = list(BH_corr)
  
  N_benef[i] = length(BH_corr[BH_corr <= FDR_critVal])
  N_tot[i] = length(p_vec)
  
  # Convert the character string to double (phenotypic effect size difference between focal genotype and its muational neighbors)
  k_vec = eval(parse(text = k_list[[i]]))
  
  if (N_benef[i] == 0) {
    # NA, if the genotype has no neighboring variant with higher performance.
    kb_list[i] = NA
    kb_mean[i] = NA
    kb_std[i] = NA
  } else {
    # Else, record the differences in the median brightness measurements between the focal genotype and the neighboring variant(s)
    kb_vec = k_vec[which(BH_corr <= FDR_critVal)]
    
    kb_list[i] = list(kb_vec)
    kb_mean[i] = mean(kb_vec)
    kb_std[i] = sqrt(var(kb_vec))
  }
  
  if (is.null(p_vec) == TRUE) {
    # If the focal genotype has no neighbor variant, set bFrac as NA.
    bFrac[i] = NA
  } else {
    # bFrac: Divide the number of neighboring variant(s) with higher performance by all neighboring variant(s)
    bFrac[i] = N_benef[i]/N_tot[[i]]
  }
}

df$qval = p_adjusted
df$bFrac = bFrac
df$N_benef = N_benef
df$Neighbor = N_tot

df$kb_val = kb_list
df$kb_mean = kb_mean
df$kb_std = kb_std

### Check Values of Std - std of WT is similar to the average std calculated below
std_WT = df$std[1] # 0.1115609

### Now bin the trait performance category
ave_std = mean(df$std[! is.na(df$std)]) # 0.09325478
wt_vec = get("0")
wt_brightness = wt_vec$medianBrightness
max_brightness = max(df$medianBrightness)
min_brightness = min(df$medianBrightness)

bin_seq = seq(from = min_brightness, to = max_brightness, by = 2*bin_parameter*ave_std)

binnedArray = array()

Performance_Lv_bin = c()

bFrac_Mean = c()
bFrac_std = c()
# bFrac_nData: the number of genotypes in the bin.
bFrac_nData = c() 

mean_kb_mean = c()
mean_kb_std = c()
kb_nData = c()

for ( i in 1:(length(bin_seq)-1) ) {
  # Subset of the data within the assigned brightness range
  Template_index = which(bin_seq[i]<df$medianBrightness & df$medianBrightness<bin_seq[i+1])
  Template_bin = df[ Template_index, ]
  Performance_Lv_bin[i] = sum(bin_seq[i],bin_seq[i+1])/2

  # Remove NA: In bFrac data, report NA when the focal genotype does not have any neighbor reported in the dataset.     
  bFrac_vec = Template_bin$bFrac
  # Take the average of bFrac of all genotypes within the brightness range and variance.
  bFrac_Mean[i] = mean(bFrac_vec[!is.na(bFrac_vec)])
  bFrac_std[i] = sqrt(var(bFrac_vec[!is.na(bFrac_vec)]))
  bFrac_nData[i] = length(bFrac_vec[!is.na(bFrac_vec)])
  
  # Remove NA: In kb data, report NA when the focal genotype does not have any neighbor with higher performance reported in the dataset. 
  kb_mean_vec = Template_bin$kb_mean
  
  if (is_empty(kb_mean_vec[!is.na(kb_mean_vec)]) == FALSE) {
    mean_kb_mean[i] = mean(kb_mean_vec[!is.na(kb_mean_vec)])
    mean_kb_std[i] = sqrt(var(kb_mean_vec[!is.na(kb_mean_vec)]))
    kb_nData[i] = length(kb_mean_vec[!is.na(kb_mean_vec)])
  } else {
    mean_kb_mean[i] = 0
    mean_kb_std[i] = 0
    kb_nData[i] = 0
  }
}

binnedArray = data.frame(Performance_Lv_bin, bFrac_Mean)

binnedArray$bFrac_std = bFrac_std
binnedArray$bFrac_nData = bFrac_nData
binnedArray$bFrac_SE = bFrac_std/sqrt(bFrac_nData)

binnedArray$mean_kb_mean = mean_kb_mean
binnedArray$mean_kb_std = mean_kb_std
binnedArray$kb_nData = kb_nData

######################## Supplementary figure - un-binned result ########################
SIplt1 = ggplot( df, aes(x=medianBrightness, y=bFrac ) )
SI_nobin = SIplt1 + geom_rect(aes(xmin = 3, xmax = +Inf, ymin = -Inf, ymax = Inf), fill = "green", alpha = 0.025) + theme_bw(base_size = 10) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank() ) + 
  geom_point() + xlab("Fluorescence (log)") + ylab("Mean fraction of beneficial mutation") + theme(
    axis.title.x = element_text(size=9),
    axis.title.y = element_text(size=9)
  )

grid.arrange(SI_nobin, SI_std3, SI_std1, SI_std2, nrow = 2, ncol = 2)

######################## Supplementary figures - Full binned data ########################

SIplt2 = ggplot( binnedArray, aes(x=Performance_Lv_bin, y=bFrac_Mean, ymin = bFrac_Mean - 1*bFrac_SE, ymax = bFrac_Mean + 1*bFrac_SE ) ) # +- 1 SE
SI2 = SIplt2 + geom_rect(aes(xmin = 3, xmax = +Inf, ymin = -Inf, ymax = Inf), fill = "green", alpha = 0.025) + theme_bw(base_size = 10) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank() ) + 
  geom_pointrange() + xlab("Fluorescence (log)") + ylab("Mean fraction of beneficial mutation") + theme(
    axis.title.x = element_text(size=9),
    axis.title.y = element_text(size=9)
  )

SIplt3 = ggplot(binnedArray, aes(x=Performance_Lv_bin, y=bFrac_nData))
SI3 = SIplt3 + geom_rect(aes(xmin = 3, xmax = +Inf, ymin = -Inf, ymax = Inf), fill = "green", alpha = 0.025) +
  theme_bw(base_size = 10) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank() )+ geom_col(colour = "black", fill = "black") + 
  xlab("Fluorescence (log)") + ylab("Bin size")

grid.arrange(SI1, SI2, nrow = 2)

SIplt4 = ggplot( binnedArray, aes(x=Performance_Lv_bin, y=mean_kb_mean, ymin = mean_kb_mean - 1*mean_kb_std, ymax = mean_kb_mean + 1*mean_kb_std ) ) # +- 1 std
SI4 = SIplt4 + geom_rect(aes(xmin = 3, xmax = +Inf, ymin = -Inf, ymax = Inf), fill = "green", alpha = 0.025) + theme_bw(base_size = 10) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank() ) + 
  geom_pointrange() + xlab("Fluorescence (log)") + ylab("Mean effect size") + theme(
    axis.title.x = element_text(size=9),
    axis.title.y = element_text(size=9)
  )

SIplt5 = ggplot(binnedArray, aes(x=Performance_Lv_bin, y=kb_nData))
SI5 = SIplt5 + geom_rect(aes(xmin = 3, xmax = +Inf, ymin = -Inf, ymax = Inf), fill = "green", alpha = 0.025) +
  theme_bw(base_size = 10) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank() )+ geom_col(colour = "black", fill = "black") + 
  xlab("Fluorescence (log)") + ylab("Bin size")


grid.arrange(SI1, SI2, SI3, SI4, nrow = 2, ncol = 2)

######################## Main figures - Binned data in the functional brightness range ########################

binnedArray_Functional = binnedArray[binnedArray$Performance_Lv_bin > 3,]

plt1 = ggplot( binnedArray_Functional, aes(x=Performance_Lv_bin, y=mean_kb_mean, ymin = mean_kb_mean - 1*mean_kb_std, ymax = mean_kb_mean + 1*mean_kb_std ) ) # +- 1 std
p1 = plt1 + geom_rect(aes(xmin = 3, xmax = +Inf, ymin = -Inf, ymax = Inf), fill = "green", alpha = 0) + theme_bw(base_size = 10) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank() ) + 
  geom_pointrange() + xlab("Fluorescence (log)") + ylab("Mean effect sizes") + theme(
    axis.title.x = element_text(size=9),
    axis.title.y = element_text(size=9)
  ) +
  geom_smooth(method='lm', se=FALSE) +
  stat_regline_equation(label.x = 3.7, label.y = 0.65, aes(label = after_stat(rr.label)))

plt2 = ggplot( binnedArray_Functional, aes(x=Performance_Lv_bin, y=bFrac_Mean, ymin = bFrac_Mean - 1*bFrac_SE, ymax = bFrac_Mean + 1*bFrac_SE ) ) # +- 1 SE
p2 = plt2 + geom_rect(aes(xmin = 3, xmax = +Inf, ymin = -Inf, ymax = Inf), fill = "green", alpha = 0) + theme_bw(base_size = 10) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank() ) + 
  geom_pointrange() + xlab("Fluorescence (log)") + ylab("Mean fraction of beneficial mutation") + theme(
    axis.title.x = element_text(size=9),
    axis.title.y = element_text(size=9)
  ) +
  geom_smooth(method='lm', se=FALSE) +
  stat_regline_equation(label.x = 3.7, label.y = 0.065, aes(label = after_stat(rr.label)))

grid.arrange(p2, p1, nrow = 2)
