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

# Fly
#######
# Total D. mel sequences
nrow(fly_data)
#[1] 6629
# Ancient in real phylostratography, not error-prone
nrow(fly_data[fly_data $V3==1 & fly_data $V6 ==1 & fly_data $V7 ==1 & fly_data $V8 ==1 & fly_data $V9 ==1 & fly_data $V10 ==1 & fly_data $V11 ==1 & fly_data $V12 ==1 & fly_data $V13 ==1 & fly_data $V14 ==1 & fly_data $V15 ==1,])
#[1] 2471
# Ancient in real phylostratography although error prone
nrow(fly_data[fly_data $V3==1 & (fly_data $V6 >1 | fly_data $V7 >1 | fly_data $V8 >1 | fly_data $V9 >1 | fly_data $V10 >1 | fly_data $V11 >1 | fly_data $V12 >1 | fly_data $V13 >1 | fly_data $V14 >1 | fly_data $V15 >1),])
#[1] 318
# Young in real phylostratigraphy, error-prone
nrow(fly_data[fly_data $V3>1 & (fly_data $V6 >1 | fly_data $V7 >1 | fly_data $V8 >1 | fly_data $V9 >1 | fly_data $V10 >1 | fly_data $V11 >1 | fly_data $V12 >1 | fly_data $V13 >1 | fly_data $V14 >1 | fly_data $V15 >1),])
#[1] 1006
# Young in real phylostratography, not error-prone
nrow(fly_data[fly_data $V3>1 & fly_data $V6 ==1 & fly_data $V7 ==1 & fly_data $V8 ==1 & fly_data $V9 ==1 & fly_data $V10 ==1 & fly_data $V11 ==1 & fly_data $V12 ==1 & fly_data $V13 ==1 & fly_data $V14 ==1 & fly_data $V15 ==1,])
#[1] 2834

# Yeast
##########
#Total S. cerevisiae sequences
nrow(yeast_data)
#[1] 5878
# Ancient in real phylostratography, not error-prone
nrow(yeast_data[yeast_data$Age == 10 & yeast_data$SimAge == 10,])
#[1] 3544
# Ancient in real phylostratography although error prone
nrow(yeast_data[yeast_data$Age == 10 & yeast_data$SimAge < 10,])
#[1] 81
# Young in real phylostratigraphy, error-prone
nrow(yeast_data[yeast_data$Age < 10 & yeast_data$SimAge < 10,])
#[1] 588
# Young in real phylostratography, not error-prone
nrow(yeast_data[yeast_data$Age < 10 & yeast_data$SimAge == 10,])
#[1] 1665