####Multiorgan HMF####

library(tidyverse)
library(reshape2)

source('DynamicNetworkID/Master.R') 
source('DynamicNetworkID/m-value_functions.R')
source('DynamicNetworkID/Hartley_functions.R')
source('DynamicNetworkID/regression_functions.R')
source('DynamicNetworkID/simulation_functions.R')
source('DynamicNetworkID/scaling_functions.R')
source('DynamicNetworkID/plotting_functions.R')


###multiorgan

##read in multiorgan and set initial timepoint to 0
multiorgan.negddct <- read.table("Organ_Norm_AllData.txt", sep = "\t", header = T)
multiorgan.negddct$Age = multiorgan.negddct$Age - 8


###for select Gene-organ combinations, set all to 0 
multi.fix <- multiorgan.negddct
multi.fix[which(multi.fix$Organ == "Liver"), which(colnames(multi.fix) %in% c("Esr2","Th","Il6"))] <- 0
multi.fix[which(multi.fix$Organ == "Liver" & multi.fix$Strain == "SHR"), "Cdkn1a"] <- 0
multi.fix[which(multi.fix$Organ == "Spleen"), which(colnames(multi.fix) %in% c("Nts","Ren"))] <- 0
multi.fix[which(multi.fix$Organ == "Lung"), "Ren"] <- 0

multiorgan.negddct2 <- multi.fix


##################
####Format SHR####
##################

##filter for SHR
multiorgan.shr <- multiorgan.negddct2 %>% filter(Strain == "SHR") %>% select(-c(Smp.Gene.ID, Strain)) %>% arrange(Age) %>% arrange(Organ)
multiorgan.shr <- multiorgan.shr[,c(2,1, order(colnames(multiorgan.shr[3:ncol(multiorgan.shr)])) + 2)]

##impute and scale
multi.shr.KNN <- multiorgan.shr
multi.shr.KNN[,3:ncol(multiorgan.shr)] <- impute::impute.knn(data.matrix(multiorgan.shr[,3:ncol(multiorgan.shr)]), k=10, rowmax=0.5, colmax=0.8, maxp=1500, rng.seed=12345)$data

datHMF.shr = scale_zeroOne(multi.shr.KNN, center = "mean")
datHMF.shr[is.na(datHMF.shr)] <- 0.5  ##if scaled data was unavailable, set to midpoint of 0.5


##run HMF method for SHR
result.shr = HMF_fit(datHMF.shr)

##plot HMF errors
library(ggplot2)
ggplot(result.shr$errors, aes(x=log(lambda_value), y=log(error), color=as.factor(alpha_value))) +
  geom_point(size=5) + xlim(-16, 6) + ylim(0, 35)

##save plots of best HMF simulation
pdf("SHR-fits.pdf")
par(mfrow=c(2,3))

sim_plots(multiorgan.shr, result.shr$best_sim)
dev.off()


##write tables for errors, simulation, and parameters
write.table(result.shr$errors, "Multiorgan_SHR_Errors.txt", sep = "\t", quote = F, row.names = F)
write.table(result.shr$best_sim, "Multiorgan_SHR_Best_Sim.txt", sep = "\t", quote = F, row.names = F)
write.table(result.shr$best_params, "Multiorgan_SHR_Best_Params.txt", sep = "\t", quote = F)

##visualize parameters
pheatmap(result.shr$best_params, cluster_rows = F, cluster_cols = F)


##################
####Format WKY####
##################

##filter for WKY
multiorgan.wky <- multiorgan.negddct2 %>% filter(Strain == "WKY") %>% select(-c(Smp.Gene.ID, Strain)) %>% arrange(Age) %>% arrange(Organ)
multiorgan.wky <- multiorgan.wky[,c(2,1, order(colnames(multiorgan.wky[3:ncol(multiorgan.wky)])) + 2)]

##impute and scale
multi.wky.KNN <- multiorgan.wky
multi.wky.KNN[,3:ncol(multiorgan.wky)] <- impute::impute.knn(data.matrix(multiorgan.wky[,3:ncol(multiorgan.wky)]), k=10, rowmax=0.5, colmax=0.8, maxp=1500, rng.seed=12345)$data

datHMF.wky = scale_zeroOne(multi.wky.KNN, center = "mean")
datHMF.wky <- datHMF.wky[-which(apply(datHMF.wky,1,function(x)(sum(is.na(x)))) > 90),] ##accounts for missing data in WKY
datHMF.wky[is.na(datHMF.wky)] <- 0.5 ##if scaled data was unavailable, set to midpoint of 0.5

##run HMF method for WKY
result.wky = HMF_fit(datHMF.wky)


##plot HMF errors
library(ggplot2)
ggplot(result.wky$errors, aes(x=log(lambda_value), y=log(error), color=as.factor(alpha_value))) +
  geom_point(size=5) + xlim(-16, 6) + ylim(0, 35)


##save plots of best HMF simulation
pdf("WKY-fits.pdf")
par(mfrow=c(2,3))

sim_plots(multiorgan.wky, result.wky$best_sim)
dev.off()


##write tables for errors, simulation, and parameters
write.table(result.wky$errors, "Multiorgan_WKY_Errors.txt", sep = "\t", quote = F, row.names = F)
write.table(result.wky$best_sim, "Multiorgan_WKY_Best_Sim.txt", sep = "\t", quote = F, row.names = F)
write.table(result.wky$best_params, "Multiorgan_WKY_Best_Params.txt", sep = "\t", quote = F)


##visualize parameters
pheatmap(result.wky$best_params, cluster_rows = F, cluster_cols = F)


########################
####Add in Brainstem####
########################

##make sure colnames match up and are in the right order
BS.set <- read.table("Brainstem_Subsetted_Genes.txt", sep = "\t", header = T)
multi.set <- multiorgan.negddct2[, which(colnames(multiorgan.negddct2) %in% colnames(BS.set))]; multi.set <- multi.set[,c(1,2,3,order(colnames(multi.set[4:ncol(multi.set)])) + 3)]
BS.set <- BS.set[,which(colnames(BS.set) %in% colnames(multi.set))]

combined.set <- rbind(multi.set, BS.set)


##################
####Format SHR####
##################

combined.shr <-combined.set %>% filter(Strain == "SHR") %>% select(-Strain) %>% arrange(Age) %>% arrange(Organ)

##fix Cyp19a1, not present in RVL
combined.shr[which(combined.shr$Organ == "RVL"), "Cyp19a1"] <- 0
combined.shr <- combined.shr[,c(2,1,3:ncol(combined.shr))]

##impute and scale
combined.shr.KNN <- combined.shr
combined.shr.KNN[,3:ncol(combined.shr)] <- impute::impute.knn(data.matrix(combined.shr[,3:ncol(combined.shr)]), k=10, rowmax=0.5, colmax=0.8, maxp=1500, rng.seed=12345)$data

datHMF.comb.shr = scale_zeroOne(combined.shr.KNN, center = "mean")
datHMF.comb.shr[is.na(datHMF.comb.shr)] <- 0.5 ##if scaled data was unavailable, set to midpoint of 0.5

##run HMF method for brainstem-multiorgan SHR
result.comb.shr = HMF_fit(datHMF.comb.shr)

##plot errors
library(ggplot2)
ggplot(result.comb.shr$errors, aes(x=log(lambda_value), y=log(error), color=as.factor(alpha_value))) +
  geom_point(size=5) + xlim(-16, 6) + ylim(0, 35)

##save plots of best HMF simulation
pdf("SHR-combined-fits.pdf")
par(mfrow=c(2,3))

sim_plots(combined.shr, result.comb.shr$best_sim)
dev.off()


##write tables for errors, simulation, and parameters
write.table(result.comb.shr$errors, "Combined_SHR_Errors.txt", sep = "\t", quote = F, row.names = F)
write.table(result.comb.shr$best_sim, "Combined_SHR_Best_Sim.txt", sep = "\t", quote = F, row.names = F)
write.table(result.comb.shr$best_params, "Combined_SHR_Best_Params.txt", sep = "\t", quote = F)

#visualize best parameters
pheatmap(result.comb.shr$best_params, cluster_rows = F, cluster_cols = F)


##################
####Format WKY####
##################

combined.wky <-combined.set %>% filter(Strain == "WKY") %>% select(-Strain) %>% arrange(Age) %>% arrange(Organ)

##get rid of Cyp19a1, not present in RVL
combined.wky[which(combined.wky$Organ == "RVL"), "Cyp19a1"] <- 0
combined.wky <- combined.wky[,c(2,1,3:ncol(combined.wky))]

combined.wky.KNN <- combined.wky
combined.wky.KNN[,3:ncol(combined.wky)] <- impute::impute.knn(data.matrix(combined.wky[,3:ncol(combined.wky)]), k=10, rowmax=0.5, colmax=0.8, maxp=1500, rng.seed=12345)$data

datHMF.comb.wky = scale_zeroOne(combined.wky.KNN, center = "mean")
datHMF.comb.wky <- datHMF.comb.wky[-which(apply(datHMF.comb.wky,1,function(x)(sum(is.na(x)))) > 80),] #gets rid of the liver sample where there was nothing
datHMF.comb.wky[is.na(datHMF.comb.wky)] <- 0.5 ##if scaled data was unavailable, set to midpoint of 0.5

##run HMF method for brainstem-multiorgan WKY
result.comb.wky = HMF_fit(datHMF.comb.wky)


library(ggplot2)
ggplot(result.comb.wky$errors, aes(x=log(lambda_value), y=log(error), color=as.factor(alpha_value))) +
  geom_point(size=5) + xlim(-16, 6) + ylim(0, 35)

##save plots of best HMF simulation
pdf("WKY-combined-fits.pdf")
par(mfrow=c(2,3))

sim_plots(combined.wky, result.comb.wky$best_sim)
dev.off()


##write tables for errors, simulation, and parameters
write.table(result.comb.wky$errors, "Combined_WKY_Errors.txt", sep = "\t", quote = F, row.names = F)
write.table(result.comb.wky$best_sim, "Combined_WKY_Best_Sim.txt", sep = "\t", quote = F, row.names = F)
write.table(result.comb.wky$best_params, "Combined_WKY_Best_Params.txt", sep = "\t", quote = F)

#visualize best parameteres 
pheatmap(result.comb.wky$best_params, cluster_rows = F, cluster_cols = F)
