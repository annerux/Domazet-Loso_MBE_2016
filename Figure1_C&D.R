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

# fly data
fly_data <-read.table("~/Domazet_etal_MBE_data/dm3_mixed_map.txt",header=FALSE)
# the columns are: 
# gene_id	 prot_id	 phylostrata	 ti	 psname	 run1	 run2	 run3	 run4	 run5	 run6	 run7	 run8	 run9	 run10
# the "phylostrata" was calculated here while the run values are from Moyers and Zhang (2015)

# Note that older phylostrata have higher values in yeast_data and lower values in fly_data

# COUNTS
##########

# yeast
youngrealyeast<-nrow(yeast_data[yeast_data $Age<10,])
youngsimyeast<-nrow(yeast_data[yeast_data $SimAge<10,])
# three youngest phylostrata
youngrealyeast3<-nrow(yeast_data[yeast_data $Age<=3,])
youngsimyeast3<-nrow(yeast_data[yeast_data $SimAge<=3,])

comparisonyeast<-matrix(c(youngrealyeast, youngsimyeast, youngrealyeast3, youngsimyeast3 ), 2,2)

# fly - average over 10 simulation runs
youngrealfly<-nrow(fly_data[fly_data$V3>1,]) 
youngsimfly<-mean(sapply(c(6:15), function(x) nrow(fly_data[fly_data[,x]>1,]) ))
# three youngest phylostrata
youngrealfly3<-nrow(fly_data[fly_data$V3>=9,]) 
youngsimfly3<-mean(sapply(c(6:15), function(x) nrow(fly_data[fly_data[,x]>=9,]) ))

comparisonfly<-matrix(c(youngrealfly, youngsimfly, youngrealfly3, youngsimfly3 ), 2,2)

# PLOT
##########

barplot(comparisonyeast , beside=TRUE, col=c("grey","black"))
barplot(comparisonfly , beside=TRUE, col=c("grey","black"))

