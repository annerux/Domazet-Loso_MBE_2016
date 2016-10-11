# Oct 11 2016
# Anne-Ruxandra Carvunis

# DATA
#########

# yeast data
yeast_data<-read.table("~/Domazet_etal_MBE_data/GeneData6.txt",header=TRUE) 
yeast_data <-yeast_data[yeast_data $Age>0,]
# this data file was provided by B. Moyers 
# All values are from Carvunis et al. (2012) except SimAge which is from Moyers and Zhang (2016)
# orfs with conservation level 0 are removed from subsequent analysis

# reduced data table
# excluding those found error-prone by Moyers and Zhang (2015) 
data_reduced <-yeast_data[yeast_data $SimAge ==10,]


# CALCULATE CORRELATIONS
#########################

# ORF length
cor_orig<-cor.test(yeast_data $Age, yeast_data $Length, method="kendall")
# data:  data$Age and data$Length
# z = 39.239, p-value < 2.2e-16
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
      # tau 
# 0.3861573 

cor_moyers<-cor.test(yeast_data $SimAge, yeast_data $Length, method="kendall")
# data:  data$SimAge and data$Length
# z = 31.233, p-value < 2.2e-16
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
     # tau 
# 0.325631 

cor_diff<-cor.test(data_reduced$Age, data_reduced$Length, method="kendall")
# data:  data_reduced$Age and data_reduced$Length
# z = 26.413, p-value < 2.2e-16
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
      # tau 
# 0.2806213 


# RNA abundance
cor_orig<-cor.test(yeast_data $Age, yeast_data $Expression, method="kendall")
# data:  data$Age and data$Expression
# z = 26.545, p-value < 2.2e-16
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
      # tau 
# 0.2612145 

cor_moyers<-cor.test(yeast_data $SimAge, yeast_data $Expression, method="kendall")
# data:  data$SimAge and data$Expression
# z = 24.573, p-value < 2.2e-16
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
      # tau 
# 0.2561782 

cor_diff<-cor.test(data_reduced$Age, data_reduced$Expression, method="kendall")
# data:  data_reduced$Age and data_reduced$Expression
# z = 16.26, p-value < 2.2e-16
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
      # tau 
# 0.1726835 

# Proximity to transcription factor binding sites
cor_orig<-cor.test(yeast_data $Age, yeast_data $TFBS, method="kendall")
# data:  data$Age and data$TFBS
# z = 6.3755, p-value = 1.823e-10
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
       # tau 
# 0.07678981 

cor_moyers<-cor.test(yeast_data $SimAge, yeast_data $TFBS, method="kendall")
# data:  data$SimAge and data$TFBS
# z = 4.705, p-value = 2.539e-06
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
       # tau 
# 0.06003736 

cor_diff<-cor.test(data_reduced$Age, data_reduced$TFBS, method="kendall")
# data:  data_reduced$Age and data_reduced$TFBS
# z = 5.2434, p-value = 1.576e-07
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
       # tau 
# 0.06818262 

# Codon Adaptation Index
cor_orig<-cor.test(yeast_data $Age, yeast_data $CAI, method="kendall")
# data:  data$Age and data$CAI
# z = 31.703, p-value < 2.2e-16
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
      # tau 
# 0.3118077 

cor_moyers<-cor.test(yeast_data $SimAge, yeast_data $CAI, method="kendall")
# data:  data$SimAge and data$CAI
# z = 20.435, p-value < 2.2e-16
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
      # tau 
# 0.2129264 

cor_diff<-cor.test(data_reduced$Age, data_reduced$CAI, method="kendall")
# data:  data_reduced$Age and data_reduced$CAI
# z = 23.946, p-value < 2.2e-16
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
      # tau 
# 0.2542695 

# Purifying selectiom
cor_orig<-cor.test(yeast_data $Age, yeast_data $PurSelection, method="kendall")
# data:  data$Age and data$PurSelection
# z = 25.814, p-value < 2.2e-16
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
      # tau 
# 0.3155525 

cor_moyers<-cor.test(yeast_data $SimAge, yeast_data $PurSelection, method="kendall")
# data:  data$SimAge and data$PurSelection
# z = 20.73, p-value < 2.2e-16
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
      # tau 
# 0.2682812 

cor_diff<-cor.test(data_reduced$Age, data_reduced$PurSelection, method="kendall")
# data:  data_reduced$Age and data_reduced$PurSelection
# z = 17.803, p-value < 2.2e-16
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
      # tau 
# 0.2342874 

# AUG context
cor_orig<-cor.test(yeast_data $Age, yeast_data $AUG, method="kendall")
# data:  data$Age and data$AUG
# z = 11.081, p-value < 2.2e-16
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
      # tau 
# 0.1334689 

cor_moyers<-cor.test(yeast_data $SimAge, yeast_data $AUG, method="kendall")
# data:  data$SimAge and data$AUG
# z = 9.516, p-value < 2.2e-16
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
      # tau 
# 0.1214286 

cor_diff<-cor.test(data_reduced$Age, data_reduced$AUG, method="kendall")
# data:  data_reduced$Age and data_reduced$AUG
# z = 7.1653, p-value = 7.76e-13
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
      # tau 
# 0.0931743 
