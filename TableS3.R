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

# kruskall wallis table

kw<-c("Property", "Age Group", "Original with and without error-prone", "Original versus Simulated")
# ORF Length
for (thisage in c(1:10)){
	ori<-yeast_data[yeast_data $Age==thisage,c("Name","Length")]
	ori["source"]<-rep(1,nrow(ori))
	moyers<-yeast_data[yeast_data $SimAge==thisage,c("Name","Length")]
	moyers["source"]<-rep(2,nrow(moyers))
	diff<-data_reduced[data_reduced$Age ==thisage, c("Name","Length")]
	diff["source"]<-rep(3,nrow(diff))
	thisdataframe<-rbind(ori,diff)
	ori_diff<-kruskal.test(Length ~ source, data=thisdataframe)$p.value
	thisdataframe<-rbind(ori,moyers)
	ori_moyers<-kruskal.test(Length ~ source, data=thisdataframe)$p.value
	kw <- rbind(kw, c("Length",thisage, ori_diff, ori_moyers))
}
# RNA abundance
for (thisage in c(1:10)){
	ori<-yeast_data[yeast_data $Age==thisage,c("Name","Expression")]
	ori["source"]<-rep(1,nrow(ori))
	moyers<-yeast_data[yeast_data $SimAge==thisage,c("Name","Expression")]
	moyers["source"]<-rep(2,nrow(moyers))
	diff<-data_reduced[data_reduced$Age ==thisage, c("Name","Expression")]
	diff["source"]<-rep(3,nrow(diff))
	thisdataframe<-rbind(ori,diff)
	ori_diff<-kruskal.test(Expression ~ source, data=thisdataframe)$p.value
	thisdataframe<-rbind(ori,moyers)
	ori_moyers<-kruskal.test(Expression ~ source, data=thisdataframe)$p.value
	kw <- rbind(kw, c("Expression",thisage, ori_diff, ori_moyers))
}
# Proximity to TFBS
for (thisage in c(1:10)){
	ori<-yeast_data[yeast_data $Age==thisage,c("Name","TFBS")]
	ori["source"]<-rep(1,nrow(ori))
	moyers<-yeast_data[yeast_data $SimAge==thisage,c("Name","TFBS")]
	moyers["source"]<-rep(2,nrow(moyers))
	diff<-data_reduced[data_reduced$Age ==thisage, c("Name","TFBS")]
	diff["source"]<-rep(3,nrow(diff))
	thisdataframe<-rbind(ori,diff)
	ori_diff<-kruskal.test(TFBS ~ source, data=thisdataframe)$p.value
	thisdataframe<-rbind(ori,moyers)
	ori_moyers<-kruskal.test(TFBS ~ source, data=thisdataframe)$p.value
	kw <- rbind(kw, c("TFBS",thisage, ori_diff, ori_moyers))
}
# Codon Adaptation Index
for (thisage in c(1:10)){
	ori<-yeast_data[yeast_data $Age==thisage,c("Name","CAI")]
	ori["source"]<-rep(1,nrow(ori))
	moyers<-yeast_data[yeast_data $SimAge==thisage,c("Name","CAI")]
	moyers["source"]<-rep(2,nrow(moyers))
	diff<-data_reduced[data_reduced$Age ==thisage, c("Name","CAI")]
	diff["source"]<-rep(3,nrow(diff))
	thisdataframe<-rbind(ori,diff)
	ori_diff<-kruskal.test(CAI ~ source, data=thisdataframe)$p.value
	thisdataframe<-rbind(ori,moyers)
	ori_moyers<-kruskal.test(CAI ~ source, data=thisdataframe)$p.value
	kw <- rbind(kw, c("CAI",thisage, ori_diff, ori_moyers))
}
# AUG context
for (thisage in c(1:10)){
	ori<-yeast_data[yeast_data $Age==thisage,c("Name","AUG")]
	ori["source"]<-rep(1,nrow(ori))
	moyers<-yeast_data[yeast_data $SimAge==thisage,c("Name","AUG")]
	moyers["source"]<-rep(2,nrow(moyers))
	diff<-data_reduced[data_reduced$Age ==thisage, c("Name","AUG")]
	diff["source"]<-rep(3,nrow(diff))
	thisdataframe<-rbind(ori,diff)
	ori_diff<-kruskal.test(AUG ~ source, data=thisdataframe)$p.value
	thisdataframe<-rbind(ori,moyers)
	ori_moyers<-kruskal.test(AUG ~ source, data=thisdataframe)$p.value
	kw <- rbind(kw, c("AUG",thisage, ori_diff, ori_moyers))
}
# Purifying Selection
for (thisage in c(1:10)){
	ori<-yeast_data[yeast_data $Age==thisage,c("Name","PurSelection")]
	ori["source"]<-rep(1,nrow(ori))
	moyers<-yeast_data[yeast_data $SimAge==thisage,c("Name","PurSelection")]
	moyers["source"]<-rep(2,nrow(moyers))
	diff<-data_reduced[data_reduced$Age ==thisage, c("Name","PurSelection")]
	diff["source"]<-rep(3,nrow(diff))
	thisdataframe<-rbind(ori,diff)
	ori_diff<-kruskal.test(PurSelection ~ source, data=thisdataframe)$p.value
	thisdataframe<-rbind(ori,moyers)
	ori_moyers<-kruskal.test(PurSelection ~ source, data=thisdataframe)$p.value
	kw <- rbind(kw, c("PurSelection",thisage, ori_diff, ori_moyers))
}




