"""
This file is part of BrainMolecularAtlas.

This file makes Fig8b and 9

Copyright (c) 2021-2022 Blue Brain Project/EPFL 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

#install.packages("extrafont")  # to use Arial font
library("extrafont")
#font_import()
#fonts()
loadfonts(device = "postscript") ## for postscript()
# https://fromthebottomoftheheap.net/2013/09/09/preparing-figures-for-plos-one-with-r/

library(data.table)
library(limma)
library(qvalue)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(dplyr)
library(gplots) 
library(tidyr)
library(stringi)
library(edgeR) 
library(reshape2)

set.seed(12345)

# Volcano plots Diff Expr
#Get eb.fit function from http://www.biostat.jhsph.edu/~kkammers/software/CVproteomics/eb.fit.r (http://www.biostat.jhsph.edu/~kkammers/software/CVproteomics/eb.fit.r)
source("http://www.biostat.jhsph.edu/~kkammers/software/CVproteomics/eb.fit.r")

dt <- fread("../../data/df4r.txt")
colnames(dt)

################################################################
#################   Hasan 2020  ################################
################################################################

dt <- dt[Study == "Hasan 2020"]

colnames(dt)
unique(dt$location) #  "spinal cord" "striatum"    "hippocampus" "cerebellum"  "brainstem"   "cortex"     
unique(dt$Organism) # "mouse"
unique(dt$Age_days) # 156
unique(dt$Age_cat) # mature adult
unique(dt$condition) #  "EAE"     "control"
unique(dt$sample_id) #


#  "spinal cord" "striatum"    "hippocampus" "cerebellum"  "brainstem"   "cortex"  

##############################################################

# "spinal cord"

##############################################################


dt <- dt[location == "spinal cord"] # analogously for other locations in Fig 9

dt <- dt[,c("sample_id","gene_id_final","log_conc_uM_medNorm")]


conc_df <- dcast(melt(dt, id.vars=c("gene_id_final","sample_id")), 
                 gene_id_final~sample_id,fun.aggregate = mean)
conc_df <- data.table(conc_df)

colnames(conc_df)
nrow(conc_df)  # 7121
conc_df <-conc_df[complete.cases(conc_df), ]
nrow(conc_df)  # 7121

# check if there are any duplicated gene_id_final
conc_df$gene_id_final %>% duplicated() %>% any()  # FALSE



conc_df[,2:7] <- exp(conc_df[,2:7])
# normalize
conc_dfm <- normalizeMedianValues(as.matrix(conc_df[,2:7])) 

cd = cbind(conc_df[,1],conc_dfm )

conc_df <- cd

boxplot(cd[,2:7])

#########################################   Diff Expr An

conc_df <- data.frame(conc_df)
dimnames(conc_df)[[1]] <- conc_df[,1]
conc_df <- conc_df[,-1]

min(conc_df)

ttest <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y)
  results$p.value
}

rawpvalue = apply(conc_df, 1, ttest, grp1 = c(1:3), grp2 = c(4:6))

hist(rawpvalue)

##transform our data into log2 base.
conc_df = log2(conc_df)

#calculate the mean of each gene per control group
control = apply(conc_df[,1:3], 1, mean)

#calcuate the mean of each gene per test group
test = apply(conc_df[, 4:6], 1, mean) 

#confirming that we have a vector of numbers
class(control) 

foldchange <- test - control  # on log transformed data

hist(foldchange, xlab = "log2 Fold Change (Control vs Test)")

results = cbind(foldchange, rawpvalue)
results = as.data.frame(results)
results$probename <- rownames(results)

library(ggplot2)
#volcano = ggplot(data = results, aes(x = foldchange, y = -1*log10(rawpvalue)))
#volcano + geom_point()

results$class <- "#9bc4fd" #"background"
results[(results$rawpvalue < 0.05)&((results$foldchange >= log2(1.15)) | (results$foldchange <= log2(0.87))) ,4] <- "#435fab" #"significant"  # 1.15 and 0.87 are like in Hasan

nrow(results[results$class=='#435fab',]) 

postscript("hasan/volcano_EAEvsControl_PerLocation/spinalCord/volcano.eps", height = 6, width = 6.83,
           family = "Arial", paper = "special", onefile = FALSE,
           horizontal = FALSE)

#op <- par(mfrow=c(2,1),mar = c(5, 4, 0.05, 0.05) + 0.1)
op <- par(mfrow=c(1,1), font.lab=1, cex.lab=1.0, font.axis=1, cex.axis=1.0,
          oma = c(1,1,0,0) + 0.2,
          mar = c(1,1,0,0) + 0.2)

plot(results$foldchange, -log10(results$rawpvalue), pch=19,cex=0.5,
     xlab="log2(fold change)", ylab="-log10(p-value)", col = results$class)

abline(v=log2(1.15), col="darkgray", lty="dotted",lwd = 3)
abline(v=log2(0.87), col="darkgray", lty="dotted",lwd = 3)

abline(h=-log10(0.05), col="darkgray", lty="dotted",lwd = 3)


par(op)
dev.off()
## EPS file
embed_fonts("hasan/volcano_EAEvsControl_PerLocation/spinalCord/volcano.eps", 
            outfile = "hasan/volcano_EAEvsControl_PerLocation/spinalCord/volcano-embed.eps",options = "-dEPSCrop") # 


markers <- results[results$class=='#435fab',]
markers_withMeta <- merge(markers,conc_df,by=0)

#write.table(markers_withMeta$Row.names , file = "hasan/volcano_EAEvsControl_PerLocation/spinalCord/all_markers_withMeta.txt", quote = FALSE,sep = "\t",
#            row.names = FALSE, col.names = FALSE)


max( markers_withMeta$foldchange)
min( markers_withMeta$foldchange)



#############
# Heatmap
#############
library(ggplot2)
library(ggplotify)
library("grid")

postscript("hasan/volcano_EAEvsControl_PerLocation/spinalCord/heatmap.eps", height = 6, width = 6.83,
           family = "Arial", paper = "special", onefile = FALSE,
           horizontal = FALSE)
#op <- par(mfrow=c(2,1),mar = c(5, 4, 0.05, 0.05) + 0.1)
op <- par(mfrow=c(1,1), font.lab=1, cex.lab=1.0, font.axis=1, cex.axis=1.0,
          oma = c(1,1,0,0) + 0.2,
          mar = c(1,1,0,0) + 0.2)



#################################

library(ComplexHeatmap)
library(circlize)

markers_withMeta <- data.frame(markers_withMeta)
mat = as.matrix(markers_withMeta[,c(6:11)])
rownames(mat) = markers_withMeta$Row.names
colnames(mat) = colnames(markers_withMeta[,c(6:11)])

length(unique(markers_withMeta$Row.names))/ length(markers_withMeta$Row.names)

# drop genes with zero variance across samples
mat <- data.frame(mat)
mat <- mat[apply(mat, 1, var) != 0, ]

mat_scaled = t(apply(mat, 1, scale))
colnames(mat_scaled) <- colnames(mat)

type = c("CON","CON","CON","EAE","EAE","EAE") 
ha = HeatmapAnnotation(type = type, annotation_name_side = "left", col = list(type = c("EAE" = "orange", "CON" = "green")))

ht_list1 = Heatmap(mat_scaled, name = "Levels", row_km = 2, 
                   top_annotation = ha,  show_row_names = FALSE,
                   show_column_names = TRUE, row_title = NULL, show_row_dend = FALSE) 

#ht_list = rowAnnotation(block = anno_block(gp = gpar(fill = 2:3, col = NA)), 
#                        width = unit(2, "mm")) + ht_list1

draw(ht_list1, ht_gap = unit(5, "mm"))


par(op)
dev.off()
## EPS file
embed_fonts("hasan/volcano_EAEvsControl_PerLocation/spinalCord/heatmap.eps", 
            outfile = "hasan/volcano_EAEvsControl_PerLocation/spinalCord/heatmap-embed.eps",options = "-dEPSCrop") # 


#ht = draw(ht); row_order(ht)
#ht = draw(ht); row_dend(ht)

r.dend <- row_dend(ht_list1) 
rcl.list <- row_order(ht_list1) 
lapply(rcl.list, function(x) length(x))  #check/confirm size clusters
dim(mat)

out<-NULL
clu<-NULL

# loop to extract genes for each cluster.
for (i in 1:length(row_order(ht_list1))){
  if (i == 1) {
    clu <- t(t(row.names(mat[row_order(ht_list1)[[i]],])))
    out <- cbind(clu, paste("cluster", i, sep=""))
    colnames(out) <- c("GeneID", "Cluster")
  } else {
    clu <- t(t(row.names(mat[row_order(ht_list1)[[i]],])))
    clu <- cbind(clu, paste("cluster", i, sep=""))
    out <- rbind(out, clu)
  }
}

#write.table(out, file= "hasan/volcano_EAEvsControl_PerLocation/spinalCord/gene_clusters_spinalCord.txt", sep="\t", quote=F, row.names=FALSE)

out <- data.table(out)

unique(out$Cluster)

out_cl1 <- out[Cluster=='cluster1']
out_cl2 <- out[Cluster=='cluster2']


(length(unique(out_cl1$GeneID))+length(unique(out_cl2)))/nrow(markers_withMeta)

length(unique(out_cl1$GeneID))/length(out_cl1$GeneID)
length(unique(out_cl2$GeneID))/length(out_cl2$GeneID)


setDT(out_cl1)
setDT(out_cl2)

out_cl1[out_cl1$GeneID %in% out_cl2$GeneID]
out_cl2[out_cl2$GeneID %in% out_cl1$GeneID]



