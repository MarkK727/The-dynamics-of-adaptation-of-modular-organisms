#### This code searches for 1-mutant neighbor of each surveyed GFP variants. 

#### List of necessary packages
#install.packages("ggpubr")
#install.packages("stringr")

library(ggpubr)
library(stringr)

#### Download and import datasets from ( Sarkisyan, K. S., Bolotin, D. A., Meer, M. V., Usmanova, D. R., Mishin, A. S., Sharonov, G. V., ... & Kondrashov, F. A. (2016) )
#### Processed data sets are available at http://dx.doi.org/10.6084/m9.figshare.3102154.

df_GFP = read.table(file = 'Desktop/GFP Dataset/nucleotide_genotypes_to_brightness.tsv', sep = '\t', header = TRUE, na.string = "")
# df_GFP = read.table(file = 'nucleotide_genotypes_to_brightness.tsv', sep = '\t', header = TRUE)
df_BCtoBrightness = read.table(file = 'Desktop/GFP Dataset/barcodes_to_brightness.tsv', sep = '\t', header = TRUE, na.string = "")
# df_BCtoBrightness = read.table(file = 'barcodes_to_brightness.tsv', sep = '\t', header = TRUE)

df_GFP = as.data.frame(df_GFP)
df_BCtoBrightness = as.data.frame(df_BCtoBrightness)

#### Goal: For each row containing a mutant strain, add a column that contains the list of neighbor index ( e.g. c(2,45,68); c(1,3); c() )

# Calculate length of the list
N = length(df_GFP[[1]])

# Create a numeric vector storing the number of mutations (i.e. n, in nth-mutant)
Mut_number = rep(NA, N)
for (i in 1:N) {
  Mut_number[i] = str_count(df_GFP[i,1], ':') + 1
}

# Set Mut_number of wt as 0 (instead of NA)
Mut_number[is.na(Mut_number)] = 0

n = unique(Mut_number)
Mut_number = data.frame(Mut_number)
df_GFP[,6] = Mut_number

# Break the full dataset into sub-dataset of nth mutant classes. e.g. Dataframe "1" contains all 1-mutant from wt. "0" is wt. 
# All names (e.g. 0,1,2,... ) contained in the "Name_List"
t = 1
Name_List = c()

for (n_i in n) {
  nam <- paste("", n_i, sep = "")
  Name_List[t] = nam
  
  Neighbor_Index = which(df_GFP[,6] == n_i)
  
  dat = df_GFP[which(df_GFP[,6] == n_i),]
  assign(nam, dat)
  
  t = t+1
}

# The neighboring class of nth mutant is (n-1)th class and (n+1)th class. 
# Identify whether nth mutant class has either one of neighboring class.
l = sort(as.numeric(Name_List))
logicVector_plus = c() # Has element TRUE if nth mutant class has n+1 neighboring class in the dataset. Else, FALSE.
logicVector_minus = c() # Has element TRUE if nth mutant class has n-1 neighboring class in the dataset. Else, FALSE.

for (i in 1:length(l)) {
  logicVector_plus[i] = any(l[i]+1 == l)
  logicVector_minus[i] = any(l[i]-1 == l)
}

logicVector_minus[logicVector_minus == FALSE] = NA
logicVector_plus[logicVector_plus == FALSE] = NA

logicVector_plus = logicVector_plus*l + 1
logicVector_minus = logicVector_minus*l - 1

# This is the main loop that identifies all mutation neighbors for each mutant strain.
for (d_i in 1:length(l)) {
  
  ######## 1. We search for a neighboring mutant class of our focal genotype in nth mutant class ######## 
  
  nth = get(as.character(l[d_i]))
  nth_mProfile = as.data.frame(nth[,1])
  
  # The if statement check whether the nth mutant class has n+1 neighbor
  if (is.na(logicVector_plus[d_i]) == FALSE) {
    # nplus calls dataframe containing (n+1)th mutant class
    nplus = get(as.character(l[d_i]+1))
    # Next, analyze the mutation notation (e.g. SA301G) to identify the loci where mutation occurred. 
    nplus_mProfile = as.data.frame(nplus[,1])
    
    # k_plus is a character array containing lists of mutations from n+1 neighbors.
    k_plus = t( apply( nplus_mProfile, 1, function(x)  unlist( strsplit(x,":") ) ) )
    if (l[d_i] < 1) {
      k_plus = t( k_plus ) # For consistency of the array format (rows -> strain, column -> mutation notation)
    }
  }
  
  # The if statement check whether the nth mutant class has n-1 neighbor
  if (is.na(logicVector_minus[d_i]) == FALSE) {
    # nplus calls dataframe containing (n-1)th mutant class
    nminus = get(as.character(l[d_i]-1))
    nminus_mProfile = as.data.frame(nminus[,1])
    
    # k_minus is a character array containing lists of mutations n-1 neighbors.
    k_minus = t( apply( nminus_mProfile, 1, function(x)  unlist( strsplit(x,":") ) ) )
    if (l[d_i] <= 2) {
      k_minus = t( k_minus ) # For consistency of the array format (rows -> strain, column -> mutation notation)
    }
  }
  
  ######## 2. Conduct hypothesis testing (Mann-Whitney U-test) to get the p-values. The alternative hypothesis is "genotype_brightness < neighbor_brightness". ######## 
  
  # Create a vector containing all p-values of hypothesis tests.
  n_pval = c()
  n_k_dist = c()
  
  # In this loop, we find the the fraction of mutational neighbors of nth mutant class that has higher performance level (i.e. brightness) than the focal genotype.
  for (j in 1:nrow(nth)) {
    tot_pval = c()
    tot_k_dist = c()
    
    # The if statement check whether the nth mutant class has n+1 neighbor mutant class
    # In this loop, we find the the fraction of n+1 neighbors of nth mutant class that has higher performance level (i.e. brightness) than the focal genotype.
    if (is.na(logicVector_plus[d_i]) == FALSE) {
      # For each jth genotype in nth mutant class, identify its n+1 neighbor: If the mutant loci list (mLoci_list) of focal genotype 
      # is a subset of mutant loci list (k[i,]) of genotypes in n+1 mutant class, they are mutational neighbors of the focal genotype.    
      genotype = nth_mProfile[j,1]
      mLoci_list = unlist(strsplit(genotype[[1]],":"))
      k = k_plus
      logic_k = logical(length(k))
      
      for (index in 1:length(mLoci_list)) {
        logic_k[which(k == mLoci_list[index])] = TRUE
      }
      
      k[logic_k == FALSE] = 0
      k[logic_k == TRUE] = 1
      
      # In this loop, we find the indices of all n+1 neighbors of nth mutant class
      # If our focal genotype is the wt, all strains of 1-mutant class are its neighbor 
      if (l[d_i] == 0) {
        k[] = 1
        neighbor = which(k == 1)
        genotype_brightness = df_BCtoBrightness$brightness[ is.na(df_BCtoBrightness$nMutations) ]
      } else {
        k = apply(k, 2, as.numeric)
        neighbor = which(rowSums(k) == l[d_i]) 
        genotype_brightness = df_BCtoBrightness$brightness[ which(df_BCtoBrightness$nMutations ==genotype[[1]]) ]
      }
      
      ### Calculate the p-values
      neighbor_mProfile = nplus$nMutations[neighbor]
      pval = c()
      k_dist = c()
      
      if (length(neighbor) != 0) {
        for (i in 1:length(neighbor_mProfile)) {
          neighbor_brightness = df_BCtoBrightness$brightness[ which(df_BCtoBrightness$nMutations == neighbor_mProfile[[i]]) ]
          # Perform Mann-Whitney test: Alternative hypothesis is  "genotype_brightness < neighbor_brightness".
          result = wilcox.test(genotype_brightness, neighbor_brightness, paired = FALSE, alternative = "less")
          pval[i] = as.numeric(result[3])
          # Measure and record effect sizes of mutation on brightness
          k_dist[i] = as.numeric( median(neighbor_brightness) - median(genotype_brightness) )
        }
      }
      
      tot_pval = c(tot_pval, pval)
      tot_k_dist = c(tot_k_dist, k_dist)
      ###
    }
    
    # The if statement check whether the nth mutant class has n-1 neighbor mutant class
    # In this loop, we find the indices of all n-1 neighbors of nth mutant class
    if (is.na(logicVector_minus[d_i]) == FALSE) {
      # For each jth genotype in nth mutant class, identify its n-1 neighbor: If the mutant loci list (mLoci_list) of focal genotype 
      # is a subset of mutant loci list (k[i,]) of genotypes in n-1 mutant class, they are mutational neighbors of the focal genotype.    
      genotype = nth_mProfile[j,1]
      mLoci_list = unlist(strsplit(genotype[[1]],":"))
      k = k_minus
      logic_k = logical(length(k))
      
      for (index in 1:length(mLoci_list)) {
        logic_k[which(k == mLoci_list[index])] = TRUE
      }
      
      k[logic_k == FALSE] = 0
      k[logic_k == TRUE] = 1
      
      # In this loop, we find the indices of all n-1 neighbors of nth mutant class
      # If our focal genotype is the 1-mutant, wt is the only neighboring class 
      if (l[d_i] <= 1) {
        k[] = 1
        neighbor = which(k == 1)
      } else {
        k = apply(k, 2, as.numeric)
        neighbor = which(rowSums(k) == l[d_i]-1)
      }
      
      ### Calculate the p-values
      neighbor_mProfile = nminus$nMutations[neighbor]
      pval = c()
      k_dist = c()
      
      # If our focal genotype is the 1-mutant, wt is the only neighboring class. If our focal genotype is wt, this part of code will not run at all
      # since nminus class does not exist
      genotype_brightness = df_BCtoBrightness$brightness[ which(df_BCtoBrightness$nMutations ==genotype[[1]]) ]
      if (l[d_i] <= 1) {
        
        if (length(neighbor) != 0) {
          for (i in 1:length(neighbor_mProfile)) {
            neighbor_brightness = df_BCtoBrightness$brightness[ is.na(df_BCtoBrightness$nMutations) ]
            # Perform Mann-Whitney test: Alternative hypothesis is  "genotype_brightness < neighbor_brightness".
            result = wilcox.test(genotype_brightness, neighbor_brightness, paired = FALSE, alternative = "less")
            pval[i] = as.numeric(result[3])
            # Measure and record effect sizes of mutation on brightness
            k_dist[i] = as.numeric( median(neighbor_brightness) - median(genotype_brightness) )
          }
        }
        
      } else {
        
        if (length(neighbor) != 0) {
          for (i in 1:length(neighbor_mProfile)) {
            neighbor_brightness = df_BCtoBrightness$brightness[ which(df_BCtoBrightness$nMutations ==neighbor_mProfile[[i]]) ]
            # Perform Mann-Whitney test: Alternative hypothesis is  "genotype_brightness < neighbor_brightness".
            result = wilcox.test(genotype_brightness, neighbor_brightness, paired = FALSE, alternative = "less")
            pval[i] = as.numeric(result[3])
            # Measure and record effect sizes of mutation on brightness
            k_dist[i] = as.numeric( median(neighbor_brightness) - median(genotype_brightness) )
          }
        }
      }
      
      tot_pval = c(tot_pval, pval)
      tot_k_dist = c(tot_k_dist, k_dist)
      ###
    }
    
    n_pval[j] = list(tot_pval)
    n_k_dist[j] = list(tot_k_dist)
    
  }
  #
  
  
  nth$pval = n_pval
  nth$k = n_k_dist
  
  assign(as.character(l[d_i]), nth)
  
}

total = get("0")
for (i in l) {
  if (i != 0) {
    total = rbind( total, get(as.character(i)) )
  }
}

### Save the dataframe "total" storing the main results for further analysis.
df <- apply(total,2,as.character)
write.csv(df, "Desktop/GFP Dataset/df.avGFP.csv", row.names=FALSE)
