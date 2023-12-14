#!/usr/bin/Rscript

list.of.packages <- c( "R.utils","tidyverse",'ggplot2','dplyr','readr','zeallot')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)




library(tidyverse)
library(dplyr)
library(ggplot2)
#library("vioplot")
library(R.utils)
library(readr)








#need -n_datasets= number of datasets parsed, -dataset_1 to n= file path to dataset -sample_1= c(ordered 1-Ntimepoint after methylation or siRNA used)
#use optionals for filtering
#-c1_prmt=int: change in 1 datasets of PRMT default 0.2,
#....

options <- R.utils::commandArgs(trailingOnly = TRUE,asValues=TRUE,default=
                                  list(c1_prmt=0.2,c2_prmt=0.15,padj_prmt=0.05,similarity_threshold=0.7,
                                       c1_reg= 0.2, c2_reg = 0.05, c1_ctrl=0.15,c2_ctrl=0.10))
print(options)
#setwd('/Users/martinapedna/Desktop/PHD/MRes/Rotation1/PRMT-APA/Scripts/')
#print(options$dataset_1)
#print(options['dataset_1'])


mkdirs('../Output')
mkdirs('../Output/Plots')
mkdirs('../Output/Files')

unique_featureID <- function(dataset){
  dataset <- dataset %>%  mutate(unique_featureID=paste(
    dataset$V8,
    unlist(strsplit(dataset$feature_id, split='_'))[seq(2,nrow(dataset)*3,3)],
    dataset$chr,'start',dataset$start,'end',dataset$end,
    sep='_'))
  
  return(dataset)
}

PRMT_meandif <- vector()
for (i in c(1:options$n_tp) ){
  sample_index_i <- as.character(paste0('sample_tp_',i))
  sample_info_i <- as.character(options[sample_index_i])
  PRMT_meandif <- append(PRMT_meandif,paste0('mean_diff',sample_info_i))
}



REG_meandif <- vector()
for (i in c(1:options$n_RBP) ){
  sample_index_i <- as.character(paste0('sample_RBP_',i))
  sample_info_i <- as.character(options[sample_index_i])
  REG_meandif <- append(REG_meandif,paste0('mean_diff',sample_info_i))
}



ds <- data.frame()
for (i in c(1:options$n_tp)) {
  ds_index <- as.character(paste0('ds_tp_',i))
  file_path=as.character(options[ds_index])
  
  sample_index <- as.character(paste0('sample_tp_',i))
  sample_info=as.character(options[sample_index])
  
  tmp<- read.csv(file=file_path)
  tmp$sample <- sample_info
  names_old <- colnames(tmp)
  colnames(tmp) <- c(names_old[1:10],'ctrl_1','ctrl_2','ctrl_3','condition_1','condition_2','condition_3',names_old[17:26])
  ds <- rbind(ds,tmp)
  print('prmt done')
}

for (i in c(1:options$n_RBP)) {
  ds_index <- as.character(paste0('ds_RBP_',i))
  file_path=as.character(options[ds_index])
  
  sample_index <- as.character(paste0('sample_RBP_',i))
  sample_info=as.character(options[sample_index])
  
  tmp<- read.csv(file=file_path)
  tmp$sample <- sample_info
  names_old <- colnames(tmp)
  colnames(tmp) <- c(names_old[1:10],'ctrl_1','ctrl_2','ctrl_3','condition_1','condition_2','condition_3',names_old[17:26])
  ds <- rbind(ds,tmp)
  print('sirna done')
}


ds <- unique_featureID(ds)
print(head(ds))
#a <- read.csv('/Users/martinapedna/Desktop/PHD/MRes/Rotation1/PRMT-APA/BATCH_B_SDMAi+ADMAi_timecourse_drimseq_tables/two_step_96hrs_with_mean_diff_coords.csv')


library(zeallot)

print('reshape')
print('might take a long time')
df_all <- data.frame(row.names=unique(ds$unique_featureID))
for (feature in rownames(df_all) ){
  
  names <- ds[ds$unique_featureID==feature,]$sample
  df_all[rownames(df_all)==feature,'unique_featureID'] <- unique(ds[ds$unique_featureID==feature,]$unique_featureID)
  df_all[rownames(df_all)==feature,'gene_ID'] <- unique(ds[ds$unique_featureID==feature,]$gene_id)
  df_all[rownames(df_all)==feature,'chr'] <- unique(ds[ds$unique_featureID==feature,]$chr)
  df_all[rownames(df_all)==feature,'start'] <- unique(ds[ds$unique_featureID==feature,]$start)
  df_all[rownames(df_all)==feature,'end'] <- unique(ds[ds$unique_featureID==feature,]$end)
  df_all[rownames(df_all)==feature,'strand'] <- unique(ds[ds$unique_featureID==feature,]$strand)
  
  df_all[rownames(df_all)==feature,paste0('mean_diff',names)] %<-% ds[ds$unique_featureID==feature,]$mean_diff
  df_all[rownames(df_all)==feature,paste0('ctrl_',names)] %<-% ds[ds$unique_featureID==feature,]$ctrl_1
  df_all[rownames(df_all)==feature,paste0('condition_',names)] %<-% ds[ds$unique_featureID==feature,]$condition_1
  df_all[rownames(df_all)==feature,paste0('V5_',names)] %<-% ds[ds$unique_featureID==feature,]$V5
  df_all[rownames(df_all)==feature,paste0('twostep_gene_padj_',names)] %<-%   ds[ds$unique_featureID==feature,]$twostep_gene_padj
  df_all[rownames(df_all)==feature,paste0('twostep_transcript_padj_',names)] %<-% ds[ds$unique_featureID==feature,]$twostep_transcript_padj
  df_all[rownames(df_all)==feature,paste0('lr_',names)] %<-% ds[ds$unique_featureID==feature,]$lr
  df_all[rownames(df_all)==feature,paste0('df_',names)] %<-% ds[ds$unique_featureID==feature,]$df
  df_all[rownames(df_all)==feature,paste0('pvalue_',names)] %<-% ds[ds$unique_featureID==feature,]$pvalue
  df_all[rownames(df_all)==feature,paste0('adj_pvalue_',names)] %<-% ds[ds$unique_featureID==feature,]$adj_pvalue
  
}
print('df_all')
print(head(df_all))

write_csv(df_all,'../Output/Files/Unfiltered_Formatted_All_Data_Sets_Final.csv')


df_all <- read_csv('../Outputs/Unfiltered_Formatted_All_Data_Sets_Final.csv')

print(head(df_all))
#PRMT
print('PRMT Category:')

print('change+ significance filtering')

features_to_keep <- vector()
for (i in c(1:options$n_tp) ){
  sample_index_i <- as.character(paste0('sample_tp_',i))
  sample_info_i <- as.character(options[sample_index_i])
  #ds_tmp <- df_all |> dplyr::filter(( (paste0('mean_diff',sample_info_i) > options$c1_prmt) | ( paste0('mean_diff',sample_info_i) < -options$c1_prmt )))
  ds_tmp <- dplyr::filter(df_all,(( df_all[,paste0('mean_diff',sample_info_i)] > options$c1_prmt) | ( df_all[,paste0('mean_diff',sample_info_i)] < -options$c1_prmt)  ))
  #print(ds_tmp)
  ds_tmp<- dplyr::filter(ds_tmp,ds_tmp[,paste0('twostep_transcript_padj_',sample_info_i)] <options$padj_prmt)
  #print(ds_tmp)
  for (e in c(1:options$n_tp) ){
    sample_index_e <- as.character(paste0('sample_tp_',e))
    sample_info_e <- as.character(options[sample_index_e])
    if (sample_info_i==sample_info_e) {next}
    #print('hereee')
    tmp_2 <- dplyr::filter(ds_tmp, ((
      ds_tmp[,paste0('mean_diff',sample_info_e)] >options$c2_prmt | ds_tmp[,paste0('mean_diff',sample_info_e)] < -options$c2_prmt )&
      ((ds_tmp[,paste0('mean_diff',sample_info_i)]> 0 & ds_tmp[,paste0('mean_diff',sample_info_e)] >0) |
        (ds_tmp[,paste0('mean_diff',sample_info_i)] < 0 & ds_tmp[,paste0('mean_diff',sample_info_e)]<0)) ))
    features_to_keep <- append(features_to_keep,tmp_2$unique_featureID)
  }
  
}




df_prmt_20_15 <- df_all[df_all$unique_featureID %in% features_to_keep, ] 

print('N feaures:')
nrow(df_prmt_20_15)
print('N genes:')
length(unique(df_prmt_20_15$gene_ID))


print('70% filtering:')


"features_to_keep <- vector()
for (i in c(1:options$n_RBP)) {
  
  sample_index_i <- as.character(paste0('sample_RBP_',i))
  sample_info_i <- as.character(options[sample_index_i])

  ds_tmp <- dplyr::filter(df_prmt_20_15,(is.na(df_prmt_20_15[,paste0('mean_diff',sample_info_i)])!=TRUE))
  ds_tmp <- dplyr::filter(ds_tmp, abs(ds_tmp[,'mean_diff72']) > options$similarity_threshold*(abs(ds_tmp[,paste0('mean_diff',sample_info_i)])))
  features_to_keep <- append(features_to_keep,ds_tmp$unique_featureID)
  
}
"

ds_tmp <- df_prmt_20_15

ds_tmp <- ds_tmp[rowSums(is.na(ds_tmp[,REG_meandif]))  < options$n_RBP, ]

#print(ds_tmp)
for (i in c(1:options$n_RBP)) {
  
  sample_index_i <- as.character(paste0('sample_RBP_',i))
  sample_info_i <- as.character(options[sample_index_i])
  
  #ds_tmp <- dplyr::filter(ds_tmp,(is.na(ds_tmp[,paste0('mean_diff',sample_info_i)])!=TRUE))
  ds_tmp <- dplyr::filter(ds_tmp, abs(ds_tmp[,'mean_diff72']) > options$similarity_threshold*(abs(ds_tmp[,paste0('mean_diff',sample_info_i)])))
  #features_to_keep <- append(features_to_keep,ds_tmp$unique_featureID)
  #print(ds_tmp)
}



df_prmt_20_15_70 <- ds_tmp
#df_prmt_20_15_70 <- df_prmt_20_15[df_prmt_20_15$unique_featureID %in% features_to_keep, ] 

print('N feaures:')
nrow(df_prmt_20_15_70)
print('N genes:')
length(unique(df_prmt_20_15_70$gene_ID))





df_prmt_20_15_70$RowMean <- rowMeans(df_prmt_20_15_70[,PRMT_meandif],na.rm = T)
df_prmt_20_15_70.1 <- df_prmt_20_15_70 %>% group_by(gene_ID) %>% slice_max(abs(RowMean),with_ties=FALSE)

all_other_features <- df_all[!(df_all$unique_featureID %in% df_prmt_20_15_70.1$unique_featureID) & (df_all$gene_ID %in% df_prmt_20_15_70.1$gene_ID), ]
all_other_features$RowMean <- rowMeans(all_other_features[,PRMT_meandif],na.rm = T)
print(head(all_other_features))

features_to_add <- vector()
for (gene in df_prmt_20_15_70.1$gene_ID) {
  if (df_prmt_20_15_70.1[df_prmt_20_15_70.1$gene_ID==gene,]$RowMean>0) {
    feature <- all_other_features[all_other_features$gene_ID==gene,] %>% filter(RowMean<0) %>% slice_min(RowMean,n=1,with_ties=FALSE)
    features_to_add<- append(features_to_add,feature$unique_featureID)} else {
      feature <- all_other_features[all_other_features$gene_ID==gene,] %>% filter(RowMean>0) %>%  slice_max(RowMean,n=1,with_ties=FALSE)
      features_to_add<- append(features_to_add,feature$unique_featureID)}
  
}

df_prmt_20_15_70_PO <- rbind(df_prmt_20_15_70.1,all_other_features[all_other_features$unique_featureID %in% features_to_add,])

print('N feaures final:')
nrow(df_prmt_20_15_70_PO)
print('N genes final:')
length(unique(df_prmt_20_15_70_PO$gene_ID))


write_csv(df_prmt_20_15_70_PO,'../Output/Files/PRMT_reg.csv')

print('Regulatable sites')

df_siRNA <- df_all
#c1_reg


features_to_keep <- vector()
for (i in c(1:options$n_RBP)) {
  sample_index_i <- as.character(paste0('sample_RBP_',i))
  sample_info_i <- as.character(options[sample_index_i])
  ds_tmp <- dplyr::filter(df_siRNA,is.na(df_siRNA$mean_diff72)!=TRUE)
  ds_tmp <- dplyr::filter(ds_tmp,abs(ds_tmp[,paste0('mean_diff',sample_info_i)])>options$similarity_threshold*abs(mean_diff72))
  features_to_keep <- append(features_to_keep,ds_tmp$unique_featureID)
}

df_siRNA_20 <- df_siRNA[df_siRNA$unique_featureID %in% features_to_keep, ] 


nrow(df_siRNA_20)
length(unique(df_siRNA_20$gene_ID))


df_siRNA_20_multiple <- tibble()

for (i in c(1:options$n_RBP)) {
  print('here22')
  sample_index_i <- as.character(paste0('sample_RBP_',i))
  sample_info_i <- as.character(options[sample_index_i])
  
  ds_tmp <-dplyr::filter(df_siRNA_20,(df_siRNA_20[,paste0('mean_diff',sample_info_i)]>options$c1_reg | 
                                        df_siRNA_20[,paste0('mean_diff',sample_info_i)]< -options$c1_reg)
                                        & df_siRNA_20[,paste0('twostep_transcript_padj_',sample_info_i)]<options$c2_reg)
  
  ds_tmp$tmp_var <- ds_tmp[,paste0('mean_diff',sample_info_i)]
  ds_tmp <- ds_tmp %>%  group_by(gene_ID) %>% slice_max(abs(tmp_var),n=1,with_ties=FALSE)
  ds_tmp_other <- df_all[!(df_all$unique_featureID %in% ds_tmp$unique_featureID) & (df_all$gene_ID %in% ds_tmp$gene_ID),]
  
  for (gene in ds_tmp$gene_ID) {
    ds_gene_tmp <- ds_tmp[ds_tmp$gene_ID==gene,]
    if (ds_gene_tmp[,paste0('mean_diff',sample_info_i)]>0) {
      ds_gene_tmp_other <- ds_tmp_other[ds_tmp_other$gene_ID==gene,]
      ds_gene_tmp_other <- dplyr::filter(ds_gene_tmp_other,ds_gene_tmp_other[,paste0('mean_diff',sample_info_i)]<0) 
      ds_gene_tmp_other$tmp_var <- ds_gene_tmp_other[,paste0('mean_diff',sample_info_i)]
      ds_gene_tmp_other <- ds_gene_tmp_other %>% slice_min(tmp_var,n=1,with_ties=FALSE)
      
      ds_gene_tmp_other <- dplyr::select(ds_gene_tmp_other,-tmp_var)
      ds_gene_tmp <- dplyr::select(ds_gene_tmp,-tmp_var)

      ds_gene_tmp_other$siRNA <- sample_info_i
      ds_gene_tmp$siRNA <- sample_info_i
        
      df_siRNA_20_multiple<- rbind(df_siRNA_20_multiple,ds_gene_tmp_other)
      df_siRNA_20_multiple<- rbind(df_siRNA_20_multiple,ds_gene_tmp)
      } else {
        ds_gene_tmp_other <- ds_tmp_other[ds_tmp_other$gene_ID==gene,]
        ds_gene_tmp_other <- dplyr::filter(ds_gene_tmp_other,ds_gene_tmp_other[,paste0('mean_diff',sample_info_i)]>0) 
        ds_gene_tmp_other$tmp_var <- ds_gene_tmp_other[,paste0('mean_diff',sample_info_i)]
        ds_gene_tmp_other <- ds_gene_tmp_other %>% slice_max(tmp_var,n=1,with_ties=FALSE)
        
        ds_gene_tmp_other <- dplyr::select(ds_gene_tmp_other,-tmp_var)
        ds_gene_tmp <- dplyr::select(ds_gene_tmp,-tmp_var)
        
        ds_gene_tmp_other$siRNA <- sample_info_i
        ds_gene_tmp$siRNA <- sample_info_i
        
        df_siRNA_20_multiple<- rbind(df_siRNA_20_multiple,ds_gene_tmp_other)
        df_siRNA_20_multiple<- rbind(df_siRNA_20_multiple,ds_gene_tmp)
        }
}

}

print('here')
nrow(df_siRNA_20_multiple)
length(unique(df_siRNA_20_multiple$gene_ID))
#stable(df_siRNA_20_multiple %>% group_by(gene_ID) %>% tally() >2)

Regulatable_sites_noPRMT <- df_siRNA_20_multiple[!(df_siRNA_20_multiple$gene_ID %in% df_prmt_20_15_70_PO$gene_ID),]
#Regulatable_sites_noPRMT$RowMean <- 0
Regulatable_sites_single <- tibble()
for (gene in unique(Regulatable_sites_noPRMT$gene_ID)) {
  ds_tmp <- Regulatable_sites_noPRMT[Regulatable_sites_noPRMT$gene_ID==gene,]
  ds_tmp$RowMean <- 0
  if (nrow(ds_tmp)==2) {
    #a <- ds_tmp[ds_tmp$gene_ID==gene,]
    sirnainfo <- ds_tmp[ds_tmp$gene_ID==gene,]$siRNA[1]
    #print(sirnainfo)
    #print('sirnainfo')
    #a <- paste0('mean_diff',sirnainfo[1])
    #ds_tmp <- dplyr::mutate(ds_tmp,RowMean=ifelse((gene_ID==gene & siRNA==sirnainfo[1]),ds_tmp[,a],''))
    #print(" ds_tmp[1,paste0('mean_diff',sirnainfo[1])]")
    #print( as.double(ds_tmp[1,paste0('mean_diff',sirnainfo[1])]))
    ds_tmp[1,]$RowMean <-  as.double(ds_tmp[1,paste0('mean_diff',sirnainfo[1])])
    #print('ok')
    ds_tmp[2,]$RowMean <- as.double(ds_tmp[2,paste0('mean_diff',sirnainfo[1])])
    #ds_tmp$RowMean <- ds_tmp[ds_tmp$gene_ID==gene,paste0('mean_diff',sirnainfo[1])]
    Regulatable_sites_single <- rbind(Regulatable_sites_single,ds_tmp)
    next}
  
  biggest_v <- 0
  biggest_s <- ''
  for (i in c(1:options$n_RBP)) {
    sample_index_i <- as.character(paste0('sample_RBP_',i))
    sample_info_i <- as.character(options[sample_index_i])
    val_tmp <- ds_tmp[ds_tmp$siRNA==sample_info_i,paste0('mean_diff',sample_info_i)]
    val_tmp <- max(abs(val_tmp), na.rm=T)
    if (val_tmp>biggest_v) {
      biggest_s <- sample_info_i}}
  
  ds_tmp <- ds_tmp[ds_tmp$siRNA==biggest_s,]
  #a <- paste0('mean_diff',biggest_s)
  
  #dplyr::mutate(ds_tmp,RowMean=ifelse((gene_ID==gene & siRNA==biggest_s),ds_tmp[,a],''))
  #print('ds_tmp')
  #print(ds_tmp)
  #print("ds_tmp[1,'RowMean'] <- ds_tmp[1,paste0('mean_diff',biggest_s)]")
  #print(ds_tmp[1,'RowMean'] <- ds_tmp[1,paste0('mean_diff',biggest_s)])
  ds_tmp[1,'RowMean'] <- as.double(ds_tmp[1,paste0('mean_diff',biggest_s)])
  ds_tmp[2,'RowMean'] <- as.double(ds_tmp[2,paste0('mean_diff',biggest_s)])
    #dplyr::mutate(ds_tmp,RowMean=ds_tmp[,paste0('mean_diff',biggest_s)])
  Regulatable_sites_single <- rbind(Regulatable_sites_single,ds_tmp)
}

nrow(Regulatable_sites_single)
length(unique(Regulatable_sites_single$gene_ID))
#print(Regulatable_sites_single[,88])

write_csv(Regulatable_sites_single,'../Output/Files/RBP_reg.csv')


print('Control sites')

#c1_ctrl=0.15,c2_ctrl=0.10

PRMT_ctrl <- vector()
for (i in c(1:options$n_tp) ){
  sample_index_i <- as.character(paste0('sample_tp_',i))
  sample_info_i <- as.character(options[sample_index_i])
  PRMT_ctrl <- append(PRMT_ctrl,paste0('ctrl_',sample_info_i))
}

REG_ctrl <- vector()
for (i in c(1:options$n_RBP) ){
  sample_index_i <- as.character(paste0('sample_RBP_',i))
  sample_info_i <- as.character(options[sample_index_i])
  REG_ctrl <- append(PRMT_ctrl,paste0('ctrl_',sample_info_i))
}








Ctrl <- df_all[(!(df_all$gene_ID %in% Regulatable_sites_single$gene_ID) & !(df_all$gene_ID %in% df_prmt_20_15_70_PO$gene_ID)),]

print('nfeature filtering beg:')
print(nrow(Ctrl))
print('n genes filterinf beg:')
length(unique(Ctrl$gene_ID))

df_ctrl <-  dplyr::filter(Ctrl,rowMeans(Ctrl[,c(PRMT_ctrl,REG_ctrl)],na.rm = TRUE) > options$c2_ctrl )

for (meandiff_v in c(REG_meandif,PRMT_meandif)) {
  df_ctrl <- dplyr::filter(df_ctrl,(df_ctrl[,meandiff_v]  > - options$c1_ctrl & df_ctrl[,meandiff_v] < options$c1_ctrl))}
  

print('nfeature filtering 10-15:')
print(nrow(df_ctrl))
print('n genes filterinf 10-15:')
length(unique(df_ctrl$gene_ID))


to_keep_control <- vector()
for (gene in unique(df_ctrl$gene_ID)) {
  
  #for each gene 
  #look at how many unique features
  uniqueFG <- df_ctrl[df_ctrl$gene_ID==gene,]$unique_featureID
  #if more than 2 features 
  if (length(uniqueFG)>=2) {
    biggest <- -100
    Sec_biggest <- -100
    f_ID_biggest <- character()
    f_ID_Sec_biggest <- character()
    #for each feature
    for (feature in uniqueFG) {
      mean_tmp <- mean(na.omit(as.vector(t(df_ctrl[df_ctrl$gene_ID==gene & df_ctrl$unique_featureID==feature,][c(PRMT_ctrl,REG_ctrl)]))))
      if (is.na(mean_tmp)) {
        mean_tmp <- -Inf
      }#take mean across all experiments controls
      #if site used more than the most one untill now
      #this site becomes the most used and the other site becomes the second most used
      #if not, check if this site is more used than the second most used than we update
      #otherwise move on to the next site
      if (mean_tmp>biggest) {
        Sec_biggest<-biggest
        biggest <- mean_tmp 
        f_ID_Sec_biggest<-f_ID_biggest
        f_ID_biggest <- feature}
      else {
        if (mean_tmp>Sec_biggest) {
          Sec_biggest <- mean_tmp 
          f_ID_Sec_biggest <- feature}
      }
    }
    #after looked at all sites, take the two most used stored in the f_ID_biggest and f_ID_Sec_biggest
    to_keep_control <- append(to_keep_control, c(f_ID_Sec_biggest,f_ID_biggest))
  }
  
}


df_control_pairs <- df_ctrl[df_ctrl$unique_featureID %in% to_keep_control,]

print('N sites final')
nrow(df_control_pairs)
print('N gnes final')
length(unique(df_control_pairs$gene_ID))

write_csv(df_control_pairs,'../Output/Control_sites.csv')

print('All')

df_control_pairs$category <- 'Control_sites'
Regulatable_sites_single$category <- 'Regulatable_sites'
df_prmt_20_15_70_PO$category <- 'PRMT_regulated'

df_control_pairs$siRNA <- '/'
df_prmt_20_15_70_PO$siRNA<- '/'

df_control_pairs$RowMean <- 0
#print(colnames(df_control_pairs))
#print(colnames(Regulatable_sites_single))
#print(colnames(df_prmt_20_15_70_PO))
df_all_categories <- rbind(df_control_pairs,Regulatable_sites_single,df_prmt_20_15_70_PO)
print(df_all_categories)

table(df_all_categories$category,df_all_categories$siRNA)
print('if everything fine should not print until:')
ds <- df_all_categories


for (FID in ds[ds$siRNA=='Sam68',]$unique_featureID) {
  
  if (ds[ds$unique_featureID==FID,]$mean_diffsiSam68 != ds[ds$unique_featureID==FID,]$RowMean) {
    print(FID)
    print(ds[ds$unique_featureID==FID,]$gene_ID)
  }
}


for ( gene in ds$gene_ID) {
  tmp <- ds[ds$gene_ID==gene,]
  if ((tmp[1,]$RowMean<0 & tmp[2,]$RowMean<0)|(tmp[1,]$RowMean>0 &tmp[2,]$RowMean>0)) {
    print(gene)
  }
}

print('here')

df_all_categories$Location <- 'distal'
for (gene in df_all_categories$gene_ID) {
  if (nrow(df_all_categories[df_all_categories$gene_ID==gene,])>2) {
    df_all_categories[df_all_categories$gene_ID==gene,]$Location <- NA
  }else{
    
    ifelse ((df_all_categories[df_all_categories$gene_ID==gene,]$strand=='+') ,
            {if(df_all_categories[df_all_categories$gene_ID==gene,]$start[1] >df_all_categories[df_all_categories$gene_ID==gene,]$start[2] ) {
              df_all_categories[df_all_categories$gene_ID==gene,][2,]$Location <- 'proximal'
            }else{
              df_all_categories[df_all_categories$gene_ID==gene,][1,]$Location <- 'proximal'
            }},
            {if(df_all_categories[df_all_categories$gene_ID==gene,]$start[1] >df_all_categories[df_all_categories$gene_ID==gene,]$start[2] ) {
              df_all_categories[df_all_categories$gene_ID==gene,][1,]$Location <- 'proximal'
            }else{
              df_all_categories[df_all_categories$gene_ID==gene,][2,]$Location <- 'proximal'
            }} )}
  
}

table(complete.cases(df_all_categories$Location ))

df_all_categories <- df_all_categories %>% group_by(gene_ID) %>% mutate(Effect= 
                                            ifelse(category!='Control_sites'
                                                   ,ifelse((Location=='distal' & RowMean>0)
                                                           ,'lengthening',
                                                           ifelse((Location=='proximal' & RowMean<0)
                                                                  ,'lengthening','shortening'
                                                           )
                                                   ),
                                                   'Control')
)


df_all_categories <- df_all_categories %>% mutate(gene_ID= unlist(strsplit(gene_ID,'_'))[1])


print('if all good should not print unitill:')
for (gene in unique(df_all_categories$gene_ID)) {
  temp <- df_all_categories[df_all_categories$gene_ID==gene,]
  if (temp[1,]$Effect!=temp[2,]$Effect) {
    print(gene)
  }
}

for (gene in unique(df_all_categories$gene_ID)) {
  temp <- df_all_categories[df_all_categories$gene_ID==gene,]
  if (temp[1,]$Location==temp[2,]$Location) {
    print(gene)
  }
}

print('here')



df_all_categories <- df_all_categories %>% select(-RowMean)
write_csv(df_all_categories,'../Output/Files/All_categories_DP_LS.csv')
