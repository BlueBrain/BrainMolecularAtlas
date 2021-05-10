"""
This file is part of BrainMolecularAtlas.

This file make PCA for Fig8a

Copyright (c) 2021 Blue Brain Project/EPFL 

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


#install.packages("extrafont") # to use Arial font in figures
library("extrafont")
#font_import()
#fonts()
loadfonts(device = "postscript") ## for postscript() #  https://fromthebottomoftheheap.net/2013/09/09/preparing-figures-for-plos-one-with-r/

library(data.table)
library(limma)
library(qvalue)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(dplyr)
library(gplots)   # contains the heatmap.2 package
library(tidyr)
library(stringi)
library(edgeR) 
library(reshape2)


############################################

# with locations

######### full data ############

dt = fread("../../data/df4r.txt")
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

dt <- dt[,c("sample_id","location","gene_id_final","log_conc_uM_medNorm")]


conc_df <- dcast(melt(dt, id.vars=c("gene_id_final","sample_id","location")), 
                 gene_id_final~location+sample_id,fun.aggregate = mean)
conc_df <- data.table(conc_df)

colnames(conc_df)
nrow(conc_df)
conc_df <-conc_df[complete.cases(conc_df), ]
nrow(conc_df)

# check if there are any duplicated gene_id_final
conc_df$gene_id_final %>% duplicated() %>% any()

# normalize
conc_df[,2:37] <- exp(conc_df[,2:37])

conc_dfm <- normalizeMedianValues(as.matrix(conc_df[,2:37])) #  limma

cd = cbind(conc_df[,1],conc_dfm )

cd[,2:37] <- log(cd[,2:37])

conc_df <- cd

#####################################################################
conc_dft <- data.frame(t(conc_df[,2:37]))

dftt <- setDT(conc_dft, keep.rownames = TRUE)[] #conc_dft changes as well

locations.pca <- prcomp(conc_dft[,2:length(colnames(conc_dft))])

library(ggfortify)
#locations.pca <- prcomp(irisdf, scale. = TRUE)

conc_dft <- data.frame(conc_dft)
conc_dft["sample"] <- c()

conc_dft <- data.table(conc_dft)


conc_dft$sample <- ifelse( grepl("CON",conc_dft$rn), "Control", 
                           ifelse(grepl("EAE",conc_dft$rn), "EAE","other"              
                           ))

conc_dft$location <- ifelse(grepl("cortex", conc_dft$rn, ignore.case = T), "cortex", 
                            ifelse(grepl("hippocampus", conc_dft$rn, ignore.case = T), "hippocampus",
                                   ifelse(grepl("striatum", conc_dft$rn, ignore.case = T), "striatum",
                                          ifelse(grepl("cerebellum", conc_dft$rn, ignore.case = T), "cerebellum",
                                                 ifelse(grepl("brainstem", conc_dft$rn, ignore.case = T), "brainstem",
                                                        ifelse(grepl("spinal cord", conc_dft$rn, ignore.case = T), "spinal cord","other"
                                                               
                                                        ))))))

postscript("hasan/volcano_EAEvsControl_PerLocation/pca.eps", height = 2, width = 6.83,
           family = "Arial", paper = "special", onefile = FALSE,
           horizontal = FALSE)
op <- par(mfrow=c(1,1), font.lab=1, cex.lab=1.0, font.axis=1, cex.axis=1.0,
          oma = c(2,2,0,0) + 0.5,
          mar = c(2,2,0,0) + 0.5)

autoplot(locations.pca, data = conc_dft, colour = "location",shape="sample",size=3) + theme_bw()

par(op)
dev.off()
## EPS file
embed_fonts("hasan/volcano_EAEvsControl_PerLocation/pca.eps", 
            outfile = "hasan/volcano_EAEvsControl_PerLocation/pca-embed.eps",options = "-dEPSCrop") # 


### this one is only for the legend fonts

postscript("hasan/volcano_EAEvsControl_PerLocation/pca-leg.eps", height = 4, width = 6.83,
           family = "Arial", paper = "special", onefile = FALSE,
           horizontal = FALSE)
op <- par(mfrow=c(1,1), font.lab=1, cex.lab=1.0, font.axis=1, cex.axis=1.0,
          oma = c(2,2,0,0) + 0.5,
          mar = c(2,2,0,0) + 0.5)

autoplot(locations.pca, data = conc_dft, colour = "location",shape="sample",size=3) + theme_bw()

par(op)
dev.off()
## EPS file
embed_fonts("hasan/volcano_EAEvsControl_PerLocation/pca-leg.eps", 
            outfile = "hasan/volcano_EAEvsControl_PerLocation/pca-leg-embed.eps",options = "-dEPSCrop") # 



