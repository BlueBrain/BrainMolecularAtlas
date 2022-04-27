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


# Fully based on demo from https://github.com/amcrisan/Adjutant (MIT license)


library(data.table)
library(adjutant)
library(dplyr)
library(ggplot2)
library(tidytext) 
library(Rtsne)

set.seed(416)

queries_enz1 = read.table("queries_enzymes1.txt", sep='\t')

get_data_fun <- function(query, mindate_y, maxdate_y, output) {
  df<-processSearch(query,mindate=mindate_y,maxdate=maxdate_y)
  df['querycol']<-query
  #write.table(df,file=output, append=TRUE,quote=FALSE,sep="\t",row.names = FALSE,col.names = FALSE)
}

for(j in seq(1,nrow(queries_enz1),by=1)){
  get_data_fun(toString(queries_enz1[j,1]),mindate_y = 1990,maxdate_y=2019,output = 'output.txt')
}

###################################
library("openxlsx")
df <- read.xlsx('../data/enzConc.xlsx', sheet = 1, startRow = 1, colNames = TRUE) # I did some look through and copied txt to xlsx

df$query <- NULL
df <- df[!duplicated(df), ]
df$PMID <- as.character(df$PMID)

tidy_df<-tidyCorpus(corpus = df)
tsneObj<-runTSNE(tidy_df,check_duplicates=FALSE)
df<-inner_join(df,tsneObj$Y,by="PMID")

ggplot(df,aes(x=tsneComp1,y=tsneComp2))+
  geom_point(alpha=0.2)+
  theme_bw()

optClusters <- optimalParam(df)

df<-inner_join(df,optClusters$retItems,by="PMID") %>% 
  mutate(tsneClusterStatus = ifelse(tsneCluster == 0, "not-clustered","clustered")) 

clusterNames <- df %>%
  dplyr::group_by(tsneCluster) %>%
  dplyr::summarise(medX = median(tsneComp1),
                   medY = median(tsneComp2)) %>%
  dplyr::filter(tsneCluster != 0)

ggplot(df,aes(x=tsneComp1,y=tsneComp2,group=tsneCluster))+
  geom_point(aes(colour = tsneClusterStatus),alpha=0.2)+
  geom_label(data=clusterNames,aes(x=medX,y=medY,label=tsneCluster),size=2,colour="red")+
  stat_ellipse(aes(alpha=tsneClusterStatus))+
  scale_colour_manual(values=c("black","blue"),name="cluster status")+
  scale_alpha_manual(values=c(1,0),name="cluster status")+ # denoise
  theme_bw()

clustNames<-df %>%
  group_by(tsneCluster)%>%
  mutate(tsneClusterNames = getTopTerms(clustPMID = PMID,clustValue=tsneCluster,topNVal = 2,tidyCorpus=tidy_df)) %>%
  select(PMID,tsneClusterNames) %>%
  ungroup()

df<-inner_join(df,clustNames,by=c("PMID","tsneCluster"))

clusterNames <- df %>%
  dplyr::group_by(tsneClusterNames) %>%
  dplyr::summarise(medX = median(tsneComp1),
                   medY = median(tsneComp2)) %>%
  dplyr::filter(tsneClusterNames != "Noise")

ggplot(df,aes(x=tsneComp1,y=tsneComp2,group=tsneClusterNames))+
  geom_point(aes(colour = tsneClusterStatus),alpha=0.2)+
  stat_ellipse(aes(alpha=tsneClusterStatus))+
  geom_label(data=clusterNames,aes(x=medX,y=medY,label=tsneClusterNames),size=3,colour="darkgreen",vjust="inward",hjust="inward")+
  scale_colour_manual(values=c("black","blue"),name="cluster status")+
  scale_alpha_manual(values=c(1,0),name="cluster status")+ 
  theme_bw()

clusterNames$medY[clusterNames$tsneClusterNames=="neuron-enolas"] = clusterNames$medY[clusterNames$tsneClusterNames=="neuron-enolas"] +3
clusterNames$medY[clusterNames$tsneClusterNames=="mrna-level"] = clusterNames$medY[clusterNames$tsneClusterNames=="mrna-level"] +2

ggplot(df,aes(x=tsneComp1,y=tsneComp2,group=tsneClusterNames))+
  geom_point(aes(colour = tsneClusterStatus),alpha=0.2)+
  stat_ellipse(aes(alpha=tsneClusterStatus))+
  geom_label(data=clusterNames[clusterNames$tsneClusterNames != 'Not-Clustered',],aes(x=medX,y=medY,label=tsneClusterNames),size=2,colour="darkblue")+
  scale_colour_manual(values=c("purple","grey"),name="cluster status")+
  scale_alpha_manual(values=c(1,0),name="cluster status")+ 
  theme_bw(base_size = 10)



