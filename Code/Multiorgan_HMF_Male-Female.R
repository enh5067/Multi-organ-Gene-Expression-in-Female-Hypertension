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



###Load female multiorgan data and brainstem data, combine (equivalent to start of Multiorgan_HMF.R)

multiorgan.negddct <- read.table("Organ_Norm_AllData.txt", sep = "\t", header = T)
multiorgan.negddct$Age = multiorgan.negddct$Age - 8


###for select Gene-organ combinations, set all to 0 
multi.fix <- multiorgan.negddct
multi.fix[which(multi.fix$Organ == "Liver"), which(colnames(multi.fix) %in% c("Esr2","Th","Il6"))] <- 0
multi.fix[which(multi.fix$Organ == "Liver" & multi.fix$Strain == "SHR"), "Cdkn1a"] <- 0
multi.fix[which(multi.fix$Organ == "Spleen"), which(colnames(multi.fix) %in% c("Nts","Ren"))] <- 0
multi.fix[which(multi.fix$Organ == "Lung"), "Ren"] <- 0

multiorgan.negddct2 <- multi.fix

BS.set <- read.table("Brainstem_Subsetted_Genes.txt", sep = "\t", header = T)
multi.set <- multiorgan.negddct2[, which(colnames(multiorgan.negddct2) %in% colnames(BS.set))]; multi.set <- multi.set[,c(1,2,3,order(colnames(multi.set[4:ncol(multi.set)])) + 3)]
BS.set <- BS.set[,which(colnames(BS.set) %in% colnames(multi.set))]


combined.set <- rbind(multi.set, BS.set)


####Load Male data for comparison

male <- read.table("Male_Norm_AllData.txt", sep = "\t", header = T); male <- male[,3:ncol(male)]; male$Age <- male$Age -4
male <- male[,c(2,3,1, (order(colnames(male)[4:ncol(male)]) +3))]; colnames(male)[which(colnames(male) == "Agtr1")] <- "Agtr1a"; colnames(male)[which(colnames(male) == "Gja1")] <- "Gja"


###subset male and female data for common genes 
common.genes <- table(c(colnames(male)[4:ncol(male)], colnames(combined.set[4:ncol(combined.set)]))); common.genes <- names(common.genes)[common.genes == 2]

male.common <- male[,which(colnames(male) %in% c(common.genes, "Age","Strain","Organ"))]
female.common <- combined.set[,which(colnames(combined.set) %in% c(common.genes, "Age","Strain","Organ"))]



##the following assesses Female and Male SHR and WKY separately based on the common genes##
  #for both datasets. Pipeline mirrors that of 'Mulitorgan_HMF.R'#



#########################
####Format Female SHR####
#########################

female.shr <-female.common %>% filter(Strain == "SHR") %>% select(-Strain) %>% arrange(Age) %>% arrange(Organ)


female.shr <- female.shr[,c(2,1,3:ncol(female.shr))]

female.shr.KNN <- female.shr
female.shr.KNN[,3:ncol(female.shr)] <- impute::impute.knn(data.matrix(female.shr[,3:ncol(female.shr)]), k=10, rowmax=0.5, colmax=0.8, maxp=1500, rng.seed=12345)$data

datHMF.female.shr = scale_zeroOne1(female.shr.KNN, center = "mean")
datHMF.female.shr[is.na(datHMF.female.shr)] <- 0.5


result.female.shr = HMF_fit(datHMF.female.shr)


library(ggplot2)
ggplot(result.female.shr$errors, aes(x=log(lambda_value), y=log(error), color=as.factor(alpha_value))) +
  geom_point(size=5) + xlim(-16, 6) + ylim(0, 35)


pdf("SHR-female-fits.pdf")
par(mfrow=c(2,3))

sim_plots(female.shr, result.female.shr$best_sim)
dev.off()

write.table(result.female.shr$errors, "Female_SHR_Errors.txt", sep = "\t", quote = F, row.names = F)
write.table(result.female.shr$best_sim, "Female_SHR_Best_Sim.txt", sep = "\t", quote = F, row.names = F)
write.table(result.female.shr$best_params, "Female_SHR_Best_Params.txt", sep = "\t", quote = F)

pheatmap(result.female.shr$best_params, cluster_rows = F, cluster_cols = F)


#########################
####Format Female WKY####
#########################

female.wky <-female.common %>% filter(Strain == "WKY") %>% select(-Strain) %>% arrange(Age) %>% arrange(Organ)

female.wky <- female.wky[,c(2,1,3:ncol(female.wky))]

female.wky.KNN <- female.wky
female.wky.KNN[,3:ncol(female.wky)] <- impute::impute.knn(data.matrix(female.wky[,3:ncol(female.wky)]), k=10, rowmax=0.5, colmax=0.8, maxp=1500, rng.seed=12345)$data

datHMF.female.wky = scale_zeroOne1(female.wky.KNN, center = "mean")
datHMF.female.wky <- datHMF.female.wky[-which(apply(datHMF.female.wky,1,function(x)(sum(is.na(x)))) > 14),] #gets rid of the liver sample where there was nothing
datHMF.female.wky[is.na(datHMF.female.wky)] <- 0.5


result.female.wky = HMF_fit(datHMF.female.wky)


library(ggplot2)
ggplot(result.female.wky$errors, aes(x=log(lambda_value), y=log(error), color=as.factor(alpha_value))) +
  geom_point(size=5) + xlim(-16, 6) + ylim(0, 35)


pdf("WKY-female-fits.pdf")
par(mfrow=c(2,3))

sim_plots(female.wky, result.female.wky$best_sim)
dev.off()

write.table(result.female.wky$errors, "Female_WKY_Errors.txt", sep = "\t", quote = F, row.names = F)
write.table(result.female.wky$best_sim, "Female_WKY_Best_Sim.txt", sep = "\t", quote = F, row.names = F)
write.table(result.female.wky$best_params, "Female_WKY_Best_Params.txt", sep = "\t", quote = F)

pheatmap(result.female.wky$best_params, cluster_rows = F, cluster_cols = F)


#######################
####Format Male SHR####
#######################

male.shr <-male.common %>% filter(Strain == "SHR") %>% select(-Strain) %>% arrange(Age) %>% arrange(Organ)


male.shr <- male.shr[,c(2,1,3:ncol(male.shr))]

male.shr.KNN <- male.shr
male.shr.KNN[,3:ncol(male.shr)] <- impute::impute.knn(data.matrix(male.shr[,3:ncol(male.shr)]), k=10, rowmax=0.5, colmax=0.8, maxp=1500, rng.seed=12345)$data

datHMF.male.shr = scale_zeroOne1(male.shr.KNN, center = "mean")
datHMF.male.shr[is.na(datHMF.male.shr)] <- 0.5


result.male.shr = HMF_fit(datHMF.male.shr)


library(ggplot2)
ggplot(result.male.shr$errors, aes(x=log(lambda_value), y=log(error), color=as.factor(alpha_value))) +
  geom_point(size=5) + xlim(-16, 6) + ylim(0, 35)


pdf("SHR-male-fits.pdf")
par(mfrow=c(2,3))

sim_plots(male.shr, result.male.shr$best_sim)
dev.off()

write.table(result.male.shr$errors, "Male_SHR_Errors.txt", sep = "\t", quote = F, row.names = F)
write.table(result.male.shr$best_sim, "Male_SHR_Best_Sim.txt", sep = "\t", quote = F, row.names = F)
write.table(result.male.shr$best_params, "Male_SHR_Best_Params.txt", sep = "\t", quote = F)

pheatmap(result.male.shr$best_params, cluster_rows = F, cluster_cols = F)


#########################
####Format Male WKY####
#########################

male.wky <-male.common %>% filter(Strain == "WKY") %>% select(-Strain) %>% arrange(Age) %>% arrange(Organ)

male.wky <- male.wky[,c(2,1,3:ncol(male.wky))]

male.wky.KNN <- male.wky
male.wky.KNN[,3:ncol(male.wky)] <- impute::impute.knn(data.matrix(male.wky[,3:ncol(male.wky)]), k=10, rowmax=0.5, colmax=0.8, maxp=1500, rng.seed=12345)$data

datHMF.male.wky = scale_zeroOne1(male.wky.KNN, center = "mean")
datHMF.male.wky[is.na(datHMF.male.wky)] <- 0.5


result.male.wky = HMF_fit(datHMF.male.wky)


library(ggplot2)
ggplot(result.male.wky$errors, aes(x=log(lambda_value), y=log(error), color=as.factor(alpha_value))) +
  geom_point(size=5) + xlim(-16, 6) + ylim(0, 35)


pdf("WKY-male-fits.pdf")
par(mfrow=c(2,3))

sim_plots(male.wky, result.male.wky$best_sim)
dev.off()

write.table(result.male.wky$errors, "Male_WKY_Errors.txt", sep = "\t", quote = F, row.names = F)
write.table(result.male.wky$best_sim, "Male_WKY_Best_Sim.txt", sep = "\t", quote = F, row.names = F)
write.table(result.male.wky$best_params, "Male_WKY_Best_Params.txt", sep = "\t", quote = F)

pheatmap(result.male.wky$best_params, cluster_rows = F, cluster_cols = F)
