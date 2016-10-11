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


# PREPARE DISTRIBUTIONS MATRICES
#################################

### ORF length
# orig
mean_orig<-aggregate(yeast_data[,'Length'], by=list(yeast_data $Age), FUN=mean,na.rm=TRUE)
sd_orig<-aggregate(yeast_data[,'Length'], by=list(yeast_data $Age), FUN=sd,na.rm=TRUE)
sder_orig<-sd_orig$x/sqrt(nrow(yeast_data))

#moyers
mean_moyers<-aggregate(yeast_data[,'Length'], by= list(yeast_data $SimAge), FUN=mean,na.rm=TRUE)
sd_moyers<-aggregate(yeast_data[,'Length'], by=list(yeast_data $SimAge), FUN=sd,na.rm=TRUE)
sder_moyers<-sd_moyers$x/sqrt(nrow(yeast_data))

#diff
mean_diff<-aggregate(data_reduced[,'Length'], by=list(data_reduced$Age), FUN=mean,na.rm=TRUE)
sd_diff<-aggregate(data_reduced[,'Length'], by=list(data_reduced$Age), FUN=sd,na.rm=TRUE)
sder_diff<-sd_diff$x/sqrt(nrow(data_reduced))

#matrix
m_length<-rbind(age=seq(1,10,1), real=mean_orig$x, sim=mean_moyers$x, robust=mean_diff$x, sder_orig,sder_moyers, sder_diff)

### RNA abundance
# orig
mean_orig<-aggregate(yeast_data[,'Expression'], by=list(yeast_data $Age), FUN=mean,na.rm=TRUE)
sd_orig<-aggregate(yeast_data[,'Expression'], by=list(yeast_data $Age), FUN=sd,na.rm=TRUE)
sder_orig<-sd_orig$x/sqrt(nrow(yeast_data))

#moyers
mean_moyers<-aggregate(yeast_data[,'Expression'], by=list(yeast_data $SimAge), FUN=mean,na.rm=TRUE)
sd_moyers<-aggregate(yeast_data[,'Expression'], by=list(yeast_data $SimAge), FUN=sd,na.rm=TRUE)
sder_moyers<-sd_moyers$x/sqrt(nrow(yeast_data))

#diff
mean_diff<-aggregate(data_reduced[,'Expression'], by=list(data_reduced$Age), FUN=mean,na.rm=TRUE)
sd_diff<-aggregate(data_reduced[,'Expression'], by=list(data_reduced$Age), FUN=sd,na.rm=TRUE)
sder_diff<-sd_diff$x/sqrt(nrow(data_reduced))

# matrix
m_expression<-rbind(age=seq(1,10,1), real=mean_orig$x, sim=mean_moyers$x, robust=mean_diff$x, sder_orig,sder_moyers, sder_diff)

### Proximity to TFBS
# orig
sum_orig<-aggregate(yeast_data[,'TFBS'], by=list(yeast_data $Age), FUN=sum,na.rm=TRUE)
n_orig<-table(yeast_data $Age)
fraction_orig<-sum_orig$x/n_orig
sder_orig<-sqrt(fraction_orig * (1-fraction_orig)/n_orig)

#moyers
sum_moyers<-aggregate(yeast_data[,'TFBS'], by=list(yeast_data $SimAge), FUN=sum,na.rm=TRUE)
n_moyers<-table(yeast_data $SimAge)
fraction_moyers<-sum_moyers$x/n_moyers
sder_moyers<-sqrt(fraction_moyers * (1-fraction_moyers)/n_moyers)

#diff
sum_diff<-aggregate(data_reduced[,'TFBS'], by=list(data_reduced$Age), FUN=sum,na.rm=TRUE)
n_diff<-table(data_reduced$Age)
fraction_diff<-sum_diff $x/n_diff
sder_diff<-sqrt(fraction_diff * (1-fraction_diff)/n_diff)

# matrix
m_tfbs<-rbind(age=seq(1,10,1), real= fraction_orig, sim= fraction_moyers, robust= fraction_diff, sder_orig,sder_moyers, sder_diff)

### Purifying selection
# orig
sum_orig<-aggregate(yeast_data[,'PurSelection'], by=list(yeast_data $Age), FUN=sum,na.rm=TRUE)
n_orig<-table(yeast_data $Age)
fraction_orig<-sum_orig$x/n_orig
sder_orig<-sqrt(fraction_orig * (1-fraction_orig)/n_orig)

#moyers
sum_moyers<-aggregate(yeast_data[,'PurSelection'], by=list(yeast_data $SimAge), FUN=sum,na.rm=TRUE)
n_moyers<-table(yeast_data $SimAge)
fraction_moyers<-sum_moyers$x/n_moyers
sder_moyers<-sqrt(fraction_moyers * (1-fraction_moyers)/n_moyers)

#diff
sum_diff<-aggregate(data_reduced[,'PurSelection'], by=list(data_reduced$Age), FUN=sum,na.rm=TRUE)
n_diff<-table(data_reduced$Age)
fraction_diff<-sum_diff $x/n_diff
sder_diff<-sqrt(fraction_diff * (1-fraction_diff)/n_diff)

# matrix
m_selection<-rbind(age=seq(1,10,1), real= fraction_orig, sim= fraction_moyers, robust= fraction_diff, sder_orig,sder_moyers, sder_diff)

### AUG Context
# orig
sum_orig<-aggregate(yeast_data[,'AUG'], by=list(yeast_data $Age), FUN=sum,na.rm=TRUE)
n_orig<-table(yeast_data $Age)
fraction_orig<-sum_orig$x/n_orig
sder_orig<-sqrt(fraction_orig * (1-fraction_orig)/n_orig)

#moyers
sum_moyers<-aggregate(yeast_data[,'AUG'], by=list(yeast_data $SimAge), FUN=sum,na.rm=TRUE)
n_moyers<-table(yeast_data $SimAge)
fraction_moyers<-sum_moyers$x/n_moyers
sder_moyers<-sqrt(fraction_moyers * (1-fraction_moyers)/n_moyers)

#diff
sum_diff<-aggregate(data_reduced[,'AUG'], by=list(data_reduced$Age), FUN=sum,na.rm=TRUE)
n_diff<-table(data_reduced$Age)
fraction_diff<-sum_diff $x/n_diff
sder_diff<-sqrt(fraction_diff * (1-fraction_diff)/n_diff)

# matrix
m_aug<-rbind(age=seq(1,10,1), real= fraction_orig, sim= fraction_moyers, robust= fraction_diff, sder_orig,sder_moyers, sder_diff)

### Codon Adaptation Index
# orig
mean_orig<-aggregate(yeast_data[,'CAI'], by=list(yeast_data $Age), FUN=median,na.rm=TRUE)
sd_orig<-aggregate(yeast_data[,'CAI'], by=list(yeast_data $Age), FUN=sd,na.rm=TRUE)
sder_orig<-1.253*sd_orig$x/sqrt(nrow(data))

#moyers
mean_moyers<-aggregate(yeast_data[,'CAI'], by=list(yeast_data $SimAge), FUN=median,na.rm=TRUE)
sd_moyers<-aggregate(yeast_data[,'CAI'], by=list(yeast_data $SimAge), FUN=sd,na.rm=TRUE)
sder_moyers<-1.253*sd_moyers$x/sqrt(nrow(data))

#diff
mean_diff<-aggregate(data_reduced[,'CAI'], by=list(data_reduced$Age), FUN=median,na.rm=TRUE)
sd_diff<-aggregate(data_reduced[,'CAI'], by=list(data_reduced$Age), FUN=sd,na.rm=TRUE)
sder_diff<-1.253*sd_diff$x/sqrt(nrow(data_reduced))

# matrix
m_cai<-rbind(age=seq(1,10,1), real=mean_orig$x, sim=mean_moyers$x, robust=mean_diff$x, sder_orig,sder_moyers, sder_diff)


# PLOT
###########

# Global plotting patterns
x_title<- "Age group"
categories<-c('sim','real','robust')
colors<-c("black","grey","white")
error_x<-seq(1,37,4)
error_witdh <-0.5
par(mfrow = c(4,2))

# ORF length
y_title<-"Mean coding sequence length (nucleotides)"
yrange<-c(0,1800)
barplot(m_length[categories,], beside=TRUE, ylab=y_title, xlab=x_title, col=colors, axes=FALSE, ylim = yrange)
axis(2, at = seq(0,1800,200), labels = seq(0,1800,200))
axis(1, at = seq(2.5,40,4), labels = seq(1,10,1))
arrows(error_x +0.5, m_length['sim',] - m_length['sder_moyers',], error_x +0.5,m_length['sim',] + m_length['sder_moyers',], length=0.02, angle=90, col="gray60" , code = 3, lwd = error_witdh)
arrows(error_x +1.5, m_length['real',] - m_length['sder_orig',], error_x  +1.5,m_length['real',] + m_length['sder_orig',], length=0.02, angle=90, col="gray60" , code = 3, lwd = error_witdh)
arrows(error_x +2.5, m_length['robust',] - m_length['sder_diff',], error_x  +2.5, m_length['robust',] + m_length['sder_diff',], length=0.02, angle=90, col="gray60" , code = 3, lwd = error_witdh)

# RNA abundance
y_title<-"Mean RNA abundance (# reads per nucleotide)"
yrange<-c(0,180)
barplot(m_expression[categories,], beside=TRUE, ylab=y_title, xlab=x_title, col=colors, axes=FALSE, ylim = yrange)
axis(2, at = seq(0,180,60), labels = seq(0,180,60))
axis(1, at = seq(2.5,40,4), labels = seq(1,10,1))
arrows(error_x +0.5, m_expression['sim',] - m_expression['sder_moyers',], error_x +0.5, m_expression['sim',] + m_expression['sder_moyers',], length=0.02, angle=90, col="gray60" , code = 3, lwd = error_witdh)
arrows(error_x +1.5, m_expression['real',] - m_expression['sder_orig',], error_x  +1.5, m_expression['real',] + m_expression['sder_orig',], length=0.02, angle=90, col="gray60" , code = 3, lwd = error_witdh)
arrows(error_x +2.5, m_expression['robust',] - m_expression['sder_diff',], error_x  +2.5, m_expression['robust',] + m_expression['sder_diff',], length=0.02, angle=90, col="gray60" , code = 3, lwd = error_witdh)

# Proximity to TFBS
y_title<-"Proportion in proximity of TF binding site"
yrange<-c(0,0.35)
barplot(m_tfbs[categories,], beside=TRUE, ylab=y_title, xlab=x_title, col=colors, axes=FALSE, ylim = yrange)
axis(2, at = seq(0,0.35,0.05), labels = seq(0,0.35,0.05))
axis(1, at = seq(2.5,40,4), labels = seq(1,10,1))
arrows(error_x +0.5, m_tfbs['sim',] - m_tfbs['sder_moyers',], error_x +0.5, m_tfbs['sim',] + m_tfbs['sder_moyers',], length=0.02, angle=90, col="gray60" , code = 3, lwd = error_witdh)
arrows(error_x +1.5, m_tfbs['real',] - m_tfbs['sder_orig',], error_x  +1.5, m_tfbs['real',] + m_tfbs['sder_orig',], length=0.02, angle=90, col="gray60" , code = 3, lwd = error_witdh)
arrows(error_x +2.5, m_tfbs['robust',] - m_tfbs['sder_diff',], error_x  +2.5, m_tfbs['robust',] + m_tfbs['sder_diff',], length=0.02, angle=90, col="gray60" , code = 3, lwd = error_witdh)

# Purifying selection
y_title<-"Proportion under significant purifying selection"
yrange<-c(0,0.7)
barplot(m_selection[categories,], beside=TRUE, ylab=y_title, xlab=x_title, col=colors, axes=FALSE, ylim = yrange)
axis(2, at = seq(0,0.7,0.1), labels = seq(0,0.7,0.1))
axis(1, at = seq(2.5,40,4), labels = seq(1,10,1))
arrows(error_x +0.5, m_selection['sim',] - m_selection['sder_moyers',], error_x +0.5, m_selection['sim',] + m_selection['sder_moyers',], length=0.02, angle=90, col="gray60" , code = 3, lwd = error_witdh)
arrows(error_x +1.5, m_selection['real',] - m_selection['sder_orig',], error_x  +1.5, m_selection['real',] + m_selection['sder_orig',], length=0.02, angle=90, col="gray60" , code = 3, lwd = error_witdh)
arrows(error_x +2.5, m_selection['robust',] - m_selection['sder_diff',], error_x  +2.5, m_selection['robust',] + m_selection['sder_diff',], length=0.02, angle=90, col="gray60" , code = 3, lwd = error_witdh)

# AUG context
y_title<-"Proportion with optimal AUG context"
yrange<-c(0,0.8)
barplot(m_aug[categories,], beside=TRUE, ylab=y_title, xlab=x_title, col=colors, axes=FALSE, ylim = yrange)
axis(2, at = seq(0,0.8,0.1), labels = seq(0,0.8,0.1))
axis(1, at = seq(2.5,40,4), labels = seq(1,10,1))
arrows(error_x +0.5, m_aug['sim',] - m_aug['sder_moyers',], error_x +0.5, m_aug['sim',] + m_aug['sder_moyers',], length=0.02, angle=90, col="gray60" , code = 3, lwd = error_witdh)
arrows(error_x +1.5, m_aug['real',] - m_aug['sder_orig',], error_x  +1.5, m_aug['real',] + m_aug['sder_orig',], length=0.02, angle=90, col="gray60" , code = 3, lwd = error_witdh)
arrows(error_x +2.5, m_aug['robust',] - m_aug['sder_diff',], error_x  +2.5, m_aug['robust',] + m_aug['sder_diff',], length=0.02, angle=90, col="gray60" , code = 3, lwd = error_witdh)

# Codon Adaptation Index
y_title<-"Median codon adaptation index"
yrange<-c(0,0.2)
barplot(m_cai[categories,], beside=TRUE, ylab=y_title, xlab=x_title, col=colors, axes=FALSE, ylim = yrange)
axis(2, at = seq(0,0.2,0.02), labels = seq(0,0.2,0.02))
axis(1, at = seq(2.5,40,4), labels = seq(1,10,1))
arrows(error_x +0.5, m_cai['sim',] - m_cai['sder_moyers',], error_x +0.5, m_cai['sim',] + m_cai['sder_moyers',], length=0.02, angle=90, col="gray60" , code = 3, lwd = error_witdh)
arrows(error_x +1.5, m_cai['real',] - m_cai['sder_orig',], error_x  +1.5, m_cai['real',] + m_cai['sder_orig',], length=0.02, angle=90, col="gray60" , code = 3, lwd = error_witdh)
arrows(error_x +2.5, m_cai['robust',] - m_cai['sder_diff',], error_x  +2.5, m_cai['robust',] + m_cai['sder_diff',], length=0.02, angle=90, col="gray60" , code = 3, lwd = error_witdh)


