"""
This file is part of BrainMolecularAtlas.

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

library(RISmed)
library(data.table)
library(openxlsx)

conc <- fread("../data/df4r_healthy.txt") 
colnames(conc)

unique(conc$condition)

conc <- conc[conc$condition %in% c("young","adult", "control","WT","SORT","") ]
  
locations <- unique(conc$location)
locInterest_br <- c("cerebellum","cortex", "hippocampus","striatum","brainstem","thalamus","amygdala" )

locInterest_cells <- c("astrocytes","microglia","neurons","oligodendrocytes")

conc$V1 <- NULL

colnames(conc)
colnames(conc)[30] #log_conc_uM_medNorm
conca <- aggregate(conc[, 30], list(conc$gene_id_final,conc$location), median)

colnames(conca)
colnames(conca) <- c('Gene names','Location','Median log conc uM median norm')

##### br & cells separately
conca_i_br <- conca[conca$Location %in% locInterest_br,]
conca_i_br$log2FCconc = NA

conca_i_cells <- conca[conca$Location %in% locInterest_cells,]
conca_i_cells$log2FCconc = NA


for (i in 1:nrow(conca_i_br)){ 
  
  backgroundConc = median(conca_i_br[(conca_i_br['Gene names']==conca_i_br[i,'Gene names'])&
                                       (conca_i_br['Location']!=conca_i_br[i,'Location']),'Median log conc uM median norm'])
  
  #conca_i_br[i,'log2FCconc'] = log2(conca_i_br[i,'Median log conc uM median norm']/ backgroundConc) # log2 genewise fold enrichment "one vs all others"
  conca_i_br[i,'log2FCconc'] = conca_i_br[i,'Median log conc uM median norm'] - backgroundConc # log e genewise fold enrichment "one vs all others"
}


for (i in 1:nrow(conca_i_cells)){ 
  
  backgroundConc = median(conca_i_cells[(conca_i_cells['Gene names']==conca_i_cells[i,'Gene names'])&
                                          (conca_i_cells['Location']!=conca_i_cells[i,'Location']),'Median log conc uM median norm'])
  
  #conca_i_cells[i,'log2FCconc'] = log2(conca_i_cells[i,'Median log conc uM median norm']/ backgroundConc) # log2 genewise fold enrichment "one vs all others"
  conca_i_cells[i,'log2FCconc'] = conca_i_cells[i,'Median log conc uM median norm'] - backgroundConc # log e  genewise fold enrichment "one vs all others"
}


conca_i_br <- within(conca_i_br,  queryPubmed <- paste(`Gene names`,Location, sep=" "))
conca_i_br <- conca_i_br[nchar(conca_i_br$`Gene names`) >1, ]
conca_i_br <- conca_i_br[is.finite(conca_i_br$log2FC), ]

conca_i_cells <- within(conca_i_cells,  queryPubmed <- paste(`Gene names`,Location, sep=" "))
conca_i_cells <- conca_i_cells[nchar(conca_i_cells$`Gene names`) >1, ]
conca_i_cells <- conca_i_cells[is.finite(conca_i_cells$log2FC), ]

#write.table(unique(conca_i_br$queryPubmed), "18nov2019_queryPubmed_locInterest_un_br_new.txt", sep="\t", quote = F, col.names = F,row.names = F)
#write.table(unique(conca_i_cells$queryPubmed), "18nov2019_queryPubmed_locInterest_un_cells_new.txt", sep="\t", quote = F, col.names = F,row.names = F)
write.table(unique(conca_i_br$queryPubmed), "28july2020_queryPubmed_locInterest_un_br_new.txt", sep="\t", quote = F, col.names = F,row.names = F)
write.table(unique(conca_i_cells$queryPubmed), "28july2020_queryPubmed_locInterest_un_cells_new.txt", sep="\t", quote = F, col.names = F,row.names = F)


#queries_br = fread("18nov2019_queryPubmed_locInterest_un_br_new.txt",sep='\n',header=F)
#queries_cells = fread("18nov2019_queryPubmed_locInterest_un_cells_new.txt",sep='\n',header=F)
queries_br = fread("28july2020_queryPubmed_locInterest_un_br_new.txt",sep='\n',header=F)
queries_cells = fread("28july2020_queryPubmed_locInterest_un_cells_new.txt",sep='\n',header=F)


#qmc = quantile(conca_i$`Median log conc uM median norm`, probs = seq(0, 1, by= 0.05))
qmc_br = quantile(conca_i_br$`Median log conc uM median norm`, probs = seq(0, 1, by= 0.01))
nrow(conca_i_br[conca_i_br$`Median log conc uM median norm` > qmc_br["99%"],]) # top 99% median conc # it was 95% before
top99_br <- conca_i_br[conca_i_br$`Median log conc uM median norm` > qmc_br["99%"],]
unique(top99_br$Location)
length(unique(top99_br$`Gene names`))

qmc_cells = quantile(conca_i_cells$`Median log conc uM median norm`, probs = seq(0, 1, by= 0.01))
nrow(conca_i_cells[conca_i_cells$`Median log conc uM median norm` > qmc_cells["99%"],]) # top 99% median conc # it was 95% before
top99_cells <- conca_i_cells[conca_i_cells$`Median log conc uM median norm` > qmc_cells["99%"],]
unique(top99_cells$Location)
length(unique(top99_cells$`Gene names`))

#write.table(unique(top99_br$queryPubmed), "18nov2019_top99genes_queryPubmed_br_new.txt", sep="\t", quote = F, col.names = F,row.names = F)
#write.table(unique(top99_cells$queryPubmed), "18nov2019_top99genes_queryPubmed_cells_new.txt", sep="\t", quote = F, col.names = F,row.names = F)
write.table(unique(top99_br$queryPubmed), "28july2020_top99genes_queryPubmed_br_new.txt", sep="\t", quote = F, col.names = F,row.names = F)
write.table(unique(top99_cells$queryPubmed), "28july2020_top99genes_queryPubmed_cells_new.txt", sep="\t", quote = F, col.names = F,row.names = F)

##############################
# now query
library(RISmed)
library(data.table)

#query = fread("18nov2019_top99genes_queryPubmed_cells_new.txt",sep='\n',header=F) 
#out_count_file_name = '18nov2019_allGNhitcounts_cells_new.txt'
query = fread("28july2020_top99genes_queryPubmed_cells_new.txt",sep='\n',header=F) 
out_count_file_name = '28july2020_allGNhitcounts_cells_new.txt'

for (i in 1:nrow(query)){
  Sys.sleep(1)
  r <- EUtilsSummary(query[i,], type='esearch', db='pubmed')
  countsPapers <- QueryCount(r)
  
  
  #write.table(data.frame(query[i,],countsPapers), paste("output_new_test_prep_sub", out_count_file_name, sep = "/"), append=T, sep="\t", quote = F, col.names = F,row.names = F)
  write.table(data.frame(query[i,],countsPapers), paste("28july2020_output", out_count_file_name, sep = "/"), append=T, sep="\t", quote = F, col.names = F,row.names = F)
  
  }


#query = fread("18nov2019_top99genes_queryPubmed_br_new.txt",sep='\n',header=F) 
#out_count_file_name = '18nov2019_allGNhitcounts_br_new.txt'
query = fread("28july2020_top99genes_queryPubmed_br_new.txt",sep='\n',header=F) 
out_count_file_name = '28july2020_allGNhitcounts_br_new.txt'

for (i in 1:nrow(query)){
  Sys.sleep(1)
  r <- EUtilsSummary(query[i,], type='esearch', db='pubmed')
  countsPapers <- QueryCount(r)
  
  #write.table(data.frame(query[i,],countsPapers), paste("output_new_test_prep_sub", out_count_file_name, sep = "/"), append=T, sep="\t", quote = F, col.names = F,row.names = F)
  write.table(data.frame(query[i,],countsPapers), paste("28july2020_output", out_count_file_name, sep = "/"), append=T, sep="\t", quote = F, col.names = F,row.names = F)
}

##############################
#countsPM_br <- fread('output_new_test_prep_sub/18nov2019_allGNhitcounts_br_new.txt')
#countsPM_cells <- fread('output_new_test_prep_sub/18nov2019_allGNhitcounts_cells_new.txt')
countsPM_br <- fread('28july2020_output/28july2020_allGNhitcounts_br_new.txt')
countsPM_cells <- fread('28july2020_output/28july2020_allGNhitcounts_cells_new.txt')


countsPM <- rbind(countsPM_br,countsPM_cells)
conca_i <- rbind(top99_br,top99_cells) #rbind(conca_i_br,conca_i_cells)

library(stringr)

colnames(countsPM) <- c('queryPubmed','countsHits')
toplot = merge(conca_i,countsPM,by='queryPubmed')
toplot$countsHitsNum <-as.numeric(as.character(toplot[,6])) 

#hist(as.numeric(toplot$countsHitsNum))

toplot2 = toplot[toplot['countsHitsNum']<500,]
toplot2$log2pubmed <- NA
toplot2$log2pubmed = log2(toplot2['countsHitsNum'])
toplot2 = toplot2[order(-toplot2[,5]), ]
toplot2[,8] <- unlist(toplot2[,8])
toplot3 = toplot2[is.finite(toplot2[,8]), ]
toplot3$log2pubmedMentions <-as.numeric(as.character(unlist(toplot3[,8])))
toplot3 <- toplot3[toplot3['log2pubmedMentions']>0,]
toplot3 <-  toplot3[is.finite(toplot3[,5]), ]
unique(toplot3['Location'])
nrow(unique(toplot3['Location'])) # 11

library(ggplot2)
library(ggrepel)
library(gtable)
library(cowplot) # ggdraw

theme_set(theme_cowplot())

unique(toplot3['Location'])

##514aa3 purple

cortex_plot <- ggplot(toplot3[toplot3['Location']=='cortex',],aes(log2FCconc,log2pubmedMentions,`Gene names`)) + 
  geom_point(color=ifelse(((toplot3[toplot3['Location']=='cortex','log2FCconc'] > 2)|(toplot3[toplot3['Location']=='cortex','log2FCconc'] < -2)), "#514aa3", "grey")) +
  geom_smooth(method = "glm") + 
  geom_text_repel(aes(label=ifelse( log2FCconc> 2 ,  `Gene names`, ifelse(log2FCconc < -2, `Gene names`, "")) ),size=2,  segment.size  = 0.2,segment.color = "grey50") + 
  ggtitle('Cortex') +xlab("specificity index") + ylab("log2(mentions)") + 
  scale_x_continuous(limits=c(-7,7)) +
  scale_y_continuous(limits=c(0, 10)) +
  theme_bw(base_size = 10) + theme(plot.title = element_text(hjust = 0.5, size=12), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10))

hippocampus_plot <- ggplot(toplot3[toplot3['Location']=='hippocampus',],aes(log2FCconc,log2pubmedMentions,`Gene names`)) + 
  geom_point(color=ifelse(((toplot3[toplot3['Location']=='hippocampus','log2FCconc'] > 2)|(toplot3[toplot3['Location']=='hippocampus','log2FCconc'] < -2)), "#514aa3", "grey")) + geom_smooth(method = "glm") + 
  geom_text_repel(aes(label=ifelse( log2FCconc> 2 ,  `Gene names`, ifelse(log2FCconc < -2, `Gene names`, "")) ),size=2, segment.size  = 0.2,segment.color = "grey50") + 
  ggtitle('Hippocampus') +xlab("specificity index") + ylab(NULL) + 
  scale_x_continuous(limits=c(-7,7)) +
  scale_y_continuous(limits=c(0, 10)) +
  theme_bw(base_size = 10) + theme(plot.title = element_text(hjust = 0.5, size=12), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10))


cerebellum_plot <- ggplot(toplot3[toplot3['Location']=='cerebellum',],aes(log2FCconc,log2pubmedMentions,`Gene names`)) + 
  geom_point(color=ifelse(((toplot3[toplot3['Location']=='cerebellum','log2FCconc'] > 2)|(toplot3[toplot3['Location']=='cerebellum','log2FCconc'] < -2)), "#514aa3", "grey")) + geom_smooth(method = "glm") + 
  geom_text_repel(aes(label=ifelse( log2FCconc> 2 ,  `Gene names`, ifelse(log2FCconc < -2, `Gene names`, "")) ),size=2,  segment.size  = 0.2,segment.color = "grey50"  ) + 
  ggtitle('Cerebellum') +xlab("specificity index") + ylab("log2(mentions)") + 
  scale_x_continuous(limits=c(-7,7)) +
  scale_y_continuous(limits=c(0, 10)) +
  theme_bw(base_size = 10) + theme(plot.title = element_text(hjust = 0.5, size=12), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10))


# "amygdala" 
amygdala_plot <- ggplot(toplot3[toplot3['Location']=='amygdala',],aes(log2FCconc,log2pubmedMentions,`Gene names`)) + 
  geom_point(color=ifelse(((toplot3[toplot3['Location']=='amygdala','log2FCconc'] > 2)|(toplot3[toplot3['Location']=='amygdala','log2FCconc'] < -2)), "#514aa3", "grey")) +
  geom_smooth(method = "glm") + 
  geom_text_repel(aes(label=ifelse( log2FCconc> 2 ,  `Gene names`, ifelse(log2FCconc < -2, `Gene names`, "")) ),size=2,  segment.size  = 0.2,segment.color = "grey50") + 
  ggtitle('Amygdala') +xlab("specificity index") + ylab("log2(mentions)") + 
  scale_x_continuous(limits=c(-7,7)) +
  scale_y_continuous(limits=c(0, 10)) +
  theme_bw(base_size = 10) + theme(plot.title = element_text(hjust = 0.5, size=12), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10))


# "brainstem"
brainstem_plot <- ggplot(toplot3[toplot3['Location']=='brainstem',],aes(log2FCconc,log2pubmedMentions,`Gene names`)) + 
  geom_point(color=ifelse(((toplot3[toplot3['Location']=='brainstem','log2FCconc'] > 2)|(toplot3[toplot3['Location']=='brainstem','log2FCconc'] < -2)), "#514aa3", "grey")) +
  geom_smooth(method = "glm") + 
  geom_text_repel(aes(label=ifelse( log2FCconc> 2 ,  `Gene names`, ifelse(log2FCconc < -2, `Gene names`, "")) ),size=2,  segment.size  = 0.2,segment.color = "grey50") + 
  ggtitle('Brainstem') +xlab("specificity index") + ylab("log2(mentions)") + 
  scale_x_continuous(limits=c(-7,7)) +
  scale_y_continuous(limits=c(0, 10)) +
  theme_bw(base_size = 10) + theme(plot.title = element_text(hjust = 0.5, size=12), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10))

# "striatum"  
striatum_plot <- ggplot(toplot3[toplot3['Location']=='striatum',],aes(log2FCconc,log2pubmedMentions,`Gene names`)) + 
  geom_point(color=ifelse(((toplot3[toplot3['Location']=='striatum','log2FCconc'] > 2)|(toplot3[toplot3['Location']=='striatum','log2FCconc'] < -2)), "#514aa3", "grey")) +
  geom_smooth(method = "glm") + 
  geom_text_repel(aes(label=ifelse( log2FCconc> 2 ,  `Gene names`, ifelse(log2FCconc < -2, `Gene names`, "")) ),size=2,  segment.size  = 0.2,segment.color = "grey50") + 
  ggtitle('Striatum') +xlab("specificity index") + ylab("log2(mentions)") + 
  scale_x_continuous(limits=c(-7,7)) +
  scale_y_continuous(limits=c(0, 10)) +
  theme_bw(base_size = 10) + theme(plot.title = element_text(hjust = 0.5, size=12), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10))

# "thalamus" 
thalamus_plot <- ggplot(toplot3[toplot3['Location']=='thalamus',],aes(log2FCconc,log2pubmedMentions,`Gene names`)) + 
  geom_point(color=ifelse(((toplot3[toplot3['Location']=='thalamus','log2FCconc'] > 2)|(toplot3[toplot3['Location']=='thalamus','log2FCconc'] < -2)), "#514aa3", "grey")) +
  geom_smooth(method = "glm") + 
  geom_text_repel(aes(label=ifelse( log2FCconc> 2 ,  `Gene names`, ifelse(log2FCconc < -2, `Gene names`, "")) ),size=2,  segment.size  = 0.2,segment.color = "grey50") + 
  ggtitle('Thalamus') +xlab("specificity index") + ylab("log2(mentions)") + 
  scale_x_continuous(limits=c(-7,7)) +
  scale_y_continuous(limits=c(0, 10)) +
  theme_bw(base_size = 10) + theme(plot.title = element_text(hjust = 0.5, size=12), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10))




###
neurons_plot <- ggplot(toplot3[toplot3['Location']=='neurons',],aes(log2FCconc,log2pubmedMentions,`Gene names`)) + 
  geom_point(color=ifelse(((toplot3[toplot3['Location']=='neurons','log2FCconc'] > 2)|(toplot3[toplot3['Location']=='neurons','log2FCconc'] < -2)), "#514aa3", "grey")) + geom_smooth(method = "glm") + 
  geom_text_repel(aes(label=ifelse( log2FCconc> 2 ,  `Gene names`, ifelse(log2FCconc < -2, `Gene names`, "")) ),size=2,  segment.size  = 0.2,segment.color = "grey50"  ) + 
  ggtitle('Neurons') +xlab("specificity index") + ylab("log2(mentions)") + 
  scale_x_continuous(limits=c(-7,7)) +
  scale_y_continuous(limits=c(0, 10)) +
  theme_bw(base_size = 10) + theme(plot.title = element_text(hjust = 0.5, size=12), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10))

astrocytes_plot <- ggplot(toplot3[toplot3['Location']=='astrocytes',],aes(log2FCconc,log2pubmedMentions,`Gene names`)) + 
  geom_point(color=ifelse(((toplot3[toplot3['Location']=='astrocytes','log2FCconc'] > 2)|(toplot3[toplot3['Location']=='astrocytes','log2FCconc'] < -2)), "#514aa3", "grey")) + geom_smooth(method = "glm") + 
  geom_text_repel(aes(label=ifelse( log2FCconc> 2 ,  `Gene names`, ifelse(log2FCconc < -2, `Gene names`, "")) ),size=2,  segment.size  = 0.2,segment.color = "grey50"  ) + 
  ggtitle('Astrocytes') +xlab("specificity index") + ylab(NULL) + 
  scale_x_continuous(limits=c(-7,7)) +
  scale_y_continuous(limits=c(0, 10)) +
  theme_bw(base_size = 10) + theme(plot.title = element_text(hjust = 0.5, size=12), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10))

microglia_plot <- ggplot(toplot3[toplot3['Location']=='microglia',],aes(log2FCconc,log2pubmedMentions,`Gene names`)) + 
  geom_point(color=ifelse(((toplot3[toplot3['Location']=='microglia','log2FCconc'] > 2)|(toplot3[toplot3['Location']=='microglia','log2FCconc'] < -2)), "#514aa3", "grey")) + geom_smooth(method = "glm") + 
  geom_text_repel(aes(label=ifelse( log2FCconc> 2 ,  `Gene names`, ifelse(log2FCconc < -2, `Gene names`, "")) ),size=2,  segment.size  = 0.2,segment.color = "grey50"  ) + 
  ggtitle('Microglia') +xlab("specificity index") + ylab("log2(mentions)") + 
  scale_x_continuous(limits=c(-7,7)) +
  scale_y_continuous(limits=c(0, 10)) +
  theme_bw(base_size = 10) + theme(plot.title = element_text(hjust = 0.5, size=12), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10))

oligodendrocytes_plot <- ggplot(toplot3[toplot3['Location']=='oligodendrocytes',],aes(log2FCconc,log2pubmedMentions,`Gene names`)) + 
  geom_point(color=ifelse(((toplot3[toplot3['Location']=='oligodendrocytes','log2FCconc'] > 2)|(toplot3[toplot3['Location']=='oligodendrocytes','log2FCconc'] < -2)), "#514aa3", "grey")) + geom_smooth(method = "glm") + 
  geom_text_repel(aes(label=ifelse( log2FCconc> 2 ,  `Gene names`, ifelse(log2FCconc < -2, `Gene names`, "")) ),size=2,  segment.size  = 0.2,segment.color = "grey50"  ) + 
  ggtitle('Oligodendrocytes') +xlab("specificity index") + ylab(NULL) + 
  scale_x_continuous(limits=c(-7,7)) +
  scale_y_continuous(limits=c(0, 10)) +
  theme_bw(base_size = 10) + theme(plot.title = element_text(hjust = 0.5, size=12), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10))


top_row1 <- plot_grid(cortex_plot, thalamus_plot, ncol = 2, labels = c('a','b') )
top_row2 <- plot_grid(hippocampus_plot,amygdala_plot, ncol = 2, labels = c('c','d') )

top_row3 <- plot_grid(striatum_plot, cerebellum_plot, ncol = 2, labels = c('e','f') )
top_row4 <- plot_grid(brainstem_plot, ncol = 2, labels = c('g') )

cells_row <- plot_grid(neurons_plot, astrocytes_plot, ncol = 2, labels = c('a','b' ) )
cells_row2 <- plot_grid(microglia_plot,oligodendrocytes_plot, ncol = 2, labels = c('c','d') )

plot_grid(top_row1,top_row2,top_row3,top_row4,  nrow = 4, rel_heights = c(1,1,1,1), label_y="Log2(Mentions)")
#save to svg and eps using export 
#ggsave('28july2020_br_logFCvsPubmedGrid9_v2.pdf',  width = 18, height = 26, units = "cm", device = "pdf") # The A4 size print measures 21.0 x 29.7cm
#ggsave('99percent_15july2020_br_logFCvsPubmedGrid9_v2.pdf',  width = 18, height = 26, units = "cm", device = "pdf") # The A4 size print measures 21.0 x 29.7cm

plot_grid( cells_row,  cells_row2,nrow = 4, rel_heights = c(1,1), label_y="Log2(Mentions)")
#save to svg and eps using export 
#ggsave('28july2020_cellTypes_logFCvsPubmedGrid9_v2.pdf',  width = 18, height = 26, units = "cm", device = "pdf") # The A4 size print measures 21.0 x 29.7cm
#ggsave('99percent_15july2020_cellTypes_logFCvsPubmedGrid9_v2.pdf',  width = 18, height = 26, units = "cm", device = "pdf") # The A4 size print measures 21.0 x 29.7cm

#plot_grid(top_row1,top_row2, cells_row,  cells_row2, nrow = 4, rel_heights = c(1,1,1,1), label_y="Log2(Mentions)")

#ggsave('19nov2019_logFCvsPubmedGrid9_v2.pdf',  width = 18, height = 26, units = "cm", device = "pdf") # The A4 size print measures 21.0 x 29.7cm
#ggsave('19nov2019_logFCvsPubmedGrid9_v2.png',  width = 18, height = 26, units = "cm", device = "png") # The A4 size print measures 21.0 x 29.7cm


#
purple_neurons <- unique(toplot3$`Gene names`[(toplot3[toplot3['Location']=='neurons','log2FCconc'] > 2)|
                                                (toplot3[toplot3['Location']=='neurons','log2FCconc'] < -2)])
pur_neurons_df <- data.table(purple_neurons)
pur_neurons_df$location <- 'neurons'

purple_astrocytes <- unique(toplot3$`Gene names`[(toplot3[toplot3['Location']=='astrocytes','log2FCconc'] > 2)|
                                                         (toplot3[toplot3['Location']=='astrocytes','log2FCconc'] < -2)])
pur_astrocytes_df <- data.table(purple_astrocytes)
pur_astrocytes_df$location <- 'astrocytes'

purple_microglia <- unique(toplot3$`Gene names`[(toplot3[toplot3['Location']=='microglia','log2FCconc'] > 2)|
                                                         (toplot3[toplot3['Location']=='microglia','log2FCconc'] < -2)])
pur_microglia_df <- data.table(purple_microglia)
pur_microglia_df$location <- 'microglia'

purple_oligodendrocytes <- unique(toplot3$`Gene names`[(toplot3[toplot3['Location']=='oligodendrocytes','log2FCconc'] > 2)|
                                                         (toplot3[toplot3['Location']=='oligodendrocytes','log2FCconc'] < -2)])
pur_oligo_df <- data.table(purple_oligodendrocytes)
pur_oligo_df$location <- 'oligodendrocytes'

names(pur_neurons_df) <- c("gene","location")
names(pur_astrocytes_df)<- c("gene","location")
names(pur_microglia_df)<- c("gene","location")
names(pur_oligo_df)<- c("gene","location")

###

purple_cortex <- unique(toplot3$`Gene names`[(toplot3[toplot3['Location']=='cortex','log2FCconc'] > 2)|
                                                (toplot3[toplot3['Location']=='cortex','log2FCconc'] < -2)])
pur_cortex_df <- data.table(purple_cortex)
pur_cortex_df$location <- 'cortex'
names(pur_cortex_df) <- c("gene","location")


purple_hippocampus <- unique(toplot3$`Gene names`[(toplot3[toplot3['Location']=='hippocampus','log2FCconc'] > 2)|
                                               (toplot3[toplot3['Location']=='hippocampus','log2FCconc'] < -2)])
pur_hippocampus_df <- data.table(purple_hippocampus)
pur_hippocampus_df$location <- 'hippocampus'
names(pur_hippocampus_df) <- c("gene","location")



purple_cerebellum <- unique(toplot3$`Gene names`[(toplot3[toplot3['Location']=='cerebellum','log2FCconc'] > 2)|
                                               (toplot3[toplot3['Location']=='cerebellum','log2FCconc'] < -2)])
pur_cerebellum_df <- data.table(purple_cerebellum)
pur_cerebellum_df$location <- 'cerebellum'
names(pur_cerebellum_df) <- c("gene","location")



purple_amygdala <- unique(toplot3$`Gene names`[(toplot3[toplot3['Location']=='amygdala','log2FCconc'] > 2)|
                                               (toplot3[toplot3['Location']=='amygdala','log2FCconc'] < -2)])
pur_amygdala_df <- data.table(purple_amygdala)
pur_amygdala_df$location <- 'amygdala'
names(pur_amygdala_df) <- c("gene","location")


purple_brainstem <- unique(toplot3$`Gene names`[(toplot3[toplot3['Location']=='brainstem','log2FCconc'] > 2)|
                                               (toplot3[toplot3['Location']=='brainstem','log2FCconc'] < -2)])
pur_brainstem_df <- data.table(purple_brainstem)
pur_brainstem_df$location <- 'brainstem'
names(pur_brainstem_df) <- c("gene","location")


purple_striatum <- unique(toplot3$`Gene names`[(toplot3[toplot3['Location']=='striatum','log2FCconc'] > 2)|
                                               (toplot3[toplot3['Location']=='striatum','log2FCconc'] < -2)])
pur_striatum_df <- data.table(purple_striatum)
pur_striatum_df$location <- 'striatum'
names(pur_striatum_df) <- c("gene","location")


purple_thalamus <- unique(toplot3$`Gene names`[(toplot3[toplot3['Location']=='thalamus','log2FCconc'] > 2)|
                                               (toplot3[toplot3['Location']=='thalamus','log2FCconc'] < -2)])
pur_thalamus_df <- data.table(purple_thalamus)
pur_thalamus_df$location <- 'thalamus'
names(pur_thalamus_df) <- c("gene","location")


purpleProtsAll <- rbind(pur_neurons_df,pur_astrocytes_df,pur_microglia_df,pur_oligo_df,
                        pur_cortex_df,pur_hippocampus_df,pur_cerebellum_df,pur_amygdala_df,pur_brainstem_df,pur_striatum_df,pur_thalamus_df)

#check
unique(toplot3$`Gene names`[!((toplot3$`Gene names` %in% top99_br$`Gene names`) | (toplot3$`Gene names` %in% top99_cells$`Gene names`))])

#write.table(purpleProtsAll, "28july2020_output/final_purpleProtsAll_99percent_28july2020.txt", sep="\t", quote = F, col.names = T,row.names = F)

# For Cytoscape-STRING
# proteins at > 99% conc

#unique(top99_br$queryPubmed) are in  "28july2020_top99genes_queryPubmed_br_new.txt", 
#unique(top99_cells$queryPubmed) are in "28july2020_top99genes_queryPubmed_cells_new.txt" 

length(unique(top99_br$queryPubmed))
length(unique(top99_cells$queryPubmed))

















