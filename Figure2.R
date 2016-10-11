# Oct 11 2016
# Anne-Ruxandra Carvunis
# Note that small variations from the results published in Figure 2 are expected by virtue of the random sampling of permutations 

# DATA
##########

# fly data
fly_data <-read.table("~/Domazet_etal_MBE_data/dm3_mixed_map.txt",header=FALSE)
# the columns are: 
# gene_id	 prot_id	 phylostrata	 ti	 psname	 run1	 run2	 run3	 run4	 run5	 run6	 run7	 run8	 run9	 run10
# the "phylostrata" was calculated here while the run values are from Moyers and Zhang (2015)

# FUNCTION SIMULATING ONE TRIAL
#################################

one_trial<-function(this_order){
start_list<-fly_data$V1[fly_data $V3>1]
sim_order<-my_permutations[[this_order]]
elements_found<-NULL
result<-NULL
for (i in sim_order){
	this_sim<-as.vector(fly_data$V1[fly_data[5+i]>1])
	elements_found<-unique(c(elements_found, this_sim))
	result<-c(result,length(intersect(start_list, elements_found)))
}
return(result)
}

# SELECTING 15 PERMUTATIONS OF 10 RUNS
#################################

my_permutations<-list(sample.int(10,10,replace=FALSE))
n<-15
while(length(my_permutations)<n){
	a<-sample.int(10,10,replace=FALSE)
	equals<-0
	for( j in c(1:length(my_permutations))){
		sample<-my_permutations[[j]]
		equals<-equals+sum(a==sample)
	}
	if(equals<10){my_permutations[[length(my_permutations)+1]]<-a}
}

# SATURATION ANALYSIS - AVERAGE OF 15 PERMUTATIONS
####################################################

all_results<-sapply(c(1:n),one_trial)
sumary<-cbind(apply(all_results,1, mean), apply(all_results,1,function(x) sd(x)/sqrt(n)))

# PLOT
#########

y_title<-"Number of genes found young in real phylostratigraphy 
and susceptible to BLAST artifact in real phylostratogtaphy"
x_title<-"Number of successive simulations"
yrange<-c(0,3900)
plot(c(1:10),sumary[,1], ylab=y_title, xlab=x_title, ylim = yrange, pch=2, cex.lab=0.5) 
abline(h=3840)


