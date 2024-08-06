library(dataVisEasy)

####################
####Negddct data####
####################

negddct.full <- read.table("Organ_Norm_AllData.txt", sep = "\t", header = T)
negddct.annots <- negddct.full[,1:4]; rownames(negddct.annots) <- negddct.annots$Smp.Gene.ID; negddct.annots$Strain <- factor(negddct.annots$Strain, levels = c("WKY","SHR"))
negddct <- negddct.full %>% tibble::column_to_rownames("Smp.Gene.ID") %>% select(-c("Strain","Organ","Age")) %>% as.matrix()%>% t()
###set up annotations
initiate_params()
set_annotations(negddct.annots)
set_annot_samps(c("Organ","Strain","Age"))
organ.colors <- RColorBrewer::brewer.pal(5,"Set1"); names(organ.colors) <- unique(negddct.annots$Organ)
age.cols <- rcartocolor::carto_pal(5,"Emrld"); names(age.cols) <- unique(negddct.annots$Age)
set_annot_cols(list('Strain'=c('WKY'='blue', 'SHR'='red'), 'Organ'= c(organ.colors), 'Age'=c(age.cols)))

myHeatmap(negddct, method = "pearson", linkage = "average")

myHeatmapByAnnotation(negddct,exact = F, c("Ebf1","Ace","Agt","Agtrap","Zbtb16","Rgs","Ren","Renbp","Agtra1a","Npr1"),groupings= c("Age","Strain","Organ"), groupings.gaps = c(0,1,1),treeheight.row = 0, fontsize.row = 15,)
myHeatmapByAnnotation(negddct,exact = F, c("Tnf", "Il1b", "Ltb4r", "Smad3", "Plgf", "Ccr1", "Ccr5","Il6","Il10"),groupings= c("Age","Strain","Organ"), groupings.gaps = c(0,1,1),treeheight.row = 0, fontsize.row = 15,)
