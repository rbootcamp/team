setwd("~/capstone")
pacman::p_load(tidyverse, janitor, ggfortify, DataExplorer, pheatmap, stats, pROC)
rm(list=ls())
theme_set( theme_bw() )

#Read TGCA files ####
tcga_case_data<-read_csv("TCGA drug response_Rboot_2.csv")
tcga_coad<-read_csv("TCGA_COAD_TMM_p400_common_Rboot.csv")
tcga_read<-read_csv("TCGA_READ_TMM_p400_common_Rboot.csv")
tcga_coad_read<-merge(tcga_coad,tcga_read,by="gene_id")
rm(tcga_coad, tcga_read)

#Filtered response dataset
fcil_case_data <- tcga_case_data %>%
  select(Cancer,bcr_patient_barcode, drug_name, measure_of_response) %>%  
  mutate(value=1) %>% 
  distinct() %>% 
  pivot_wider(names_from = drug_name, values_from = value) %>%  
  filter(Fluorouracil==1) %>%
  discard(~all(is.na(.)))

rm(tcga_case_data)

#Recode patient ids from expression dataset to match response data ####
patient_id<-colnames(tcga_coad_read)[-1] %>%
  gsub("\\.", "-", .)
patient_id<-str_split_fixed(patient_id, "-", 4)[,-4]
patient_id<-paste0(patient_id[,1],"-", patient_id[,2],"-", patient_id[,3] )

#Rename Patient ID in Gene Expression Matrix
colnames(tcga_coad_read)<-c("gene_id",patient_id)

#Make vector of patient ids common between response and expression datasets
common_id <- intersect(patient_id, fcil_case_data$bcr_patient_barcode) %>%
  unique()

#Remove patients without expression data
response_data <- fcil_case_data %>%
  filter(!bcr_patient_barcode%in%
           setdiff(fcil_case_data$bcr_patient_barcode,
                   common_id))

#Transpose expression data####
#Patients are rows, genes are columns
gene_matrix <- tcga_coad_read %>%
  select(c("gene_id", patient_id)) %>%
  column_to_rownames("gene_id") %>%
  t() %>%
  data.frame() %>%
  Filter(var, .) %>%
  rownames_to_column("bcr_patient_barcode")

rm(fcil_case_data, tcga_coad_read, patient_id, common_id)

#Construct combined dataset
comb <- inner_join(response_data, gene_matrix) %>%
  select( !which(colSums(.[,sapply(.,is.numeric)])==0) %>% names() ) %>% #filter columns with no reads
  mutate(response_status = ifelse(measure_of_response %in% c("Partial Response", "Complete Response"),
                                  "Responder", "Nonresponder"), .after = measure_of_response) %>%
  mutate(response_binary = recode(response_status, #add columns for binomial response status
                                  "Responder" = 1, "Nonresponder" = 0), .after = response_status)
rm(response_data, gene_matrix)

#Make vector of gene IDs
gid <- comb %>%
  select(starts_with("ENSG")) %>%
  colnames()

#Normalize####
#NORMALISATION FUNCTION#
cal_log2fc <- function(df) {
  lg <- log2(df+1)
  med <- log2(median(df+1))
  output <- lg-med
  return(output)
}

#sapply for normalisation####

nrm <- comb
nrm[,gid] <- sapply(comb[,gid], cal_log2fc )

#Data exploration and Visualisation ####

#PCA (RAW DATA)
pc <- prcomp(comb %>% select(all_of(gid)),
             scale = TRUE)

raw_PCA<-autoplot(pc, data=comb, col="response_status")

#PCA (NORMALISED DATA)
pc_nrm <- prcomp(nrm %>% select(all_of(gid)),
                 scale = TRUE)

norm_PCA<-autoplot(pc_nrm, data=nrm, col="response_status")

#FOR LOOP FOR WILCOX-TEST ####

#Make dataframe for storing p-values
pval <- gid %>%
  data.frame(gene=.) %>% 
  add_column(p=NA) %>% 
  column_to_rownames('gene')

# Wilcox-test
for(i in gid){
  pval[i,"p"] <- wilcox.test(as.numeric(as.matrix(nrm[,i])) ~
                               comb$response_status,
                               na.action = na.omit)$p.value
}

pval<-as.data.frame(pval)

#Explore significant genes ####

#BenjaminiHochberg correction

BHpval<-data.frame(rownames(pval)[which((p.adjust(pval$p, method = "BH"))<=0.05)],p.adjust(pval$p, method = "BH")[which((p.adjust(pval$p, method = "BH"))<=0.05)])
BHgenes<-BHpval[,1]
BHgenes<-BHgenes[1:16]

#HEATMAP####

#Subset BH significant genes
mapdata <- nrm %>%
  select(bcr_patient_barcode, all_of(BHgenes)) %>%
  column_to_rownames("bcr_patient_barcode")

#Correlation Matrix
plot_correlation(mapdata %>%
                   as.matrix())

#The most highly correlated genes

for (i in 1:(ncol(mapdata)-1)){
  for (j in (i+1):(ncol(mapdata))){
    print(cor.test(mapdata[,i], mapdata[,j])$p.value)
    
 }
}


rm(i)
rm(j)

#Table to annotate by response
response_table <- nrm %>%
  select(bcr_patient_barcode, response_status) %>%
  column_to_rownames("bcr_patient_barcode")

#Heatmap
pheatmap(mapdata,
         annotation_row = response_table,
         show_rownames=F,
         fontsize=7)

# GLM ####

#Identify genes with high STERR
geneSTERR<-nrm %>%
  select(bcr_patient_barcode, response_binary, all_of(BHgenes)) %>%
  column_to_rownames("bcr_patient_barcode") %>%
  glm(response_binary~.,
      data = ., family = "binomial")
 
STERRmatrix<-summary(geneSTERR)$coef

which(STERRmatrix$Estimate-STERRmatrix$"Std. Error")


STERRmatrix %>% as.data.frame() %>%  mutate(diff= "Estimate"-"Std.Error")

#Remove Genes with high STERR

highSTERR<-c("ENSG00000105507", "ENSG00000207395", "ENSG00000230355", "ENSG00000231702", "ENSG00000250956", "ENSG00000252297", "ENSG00000253152" )

goodGENES<- setdiff(BHgenes, highSTERR)

modeldata <- nrm %>%
  select(bcr_patient_barcode, response_binary, all_of(goodGENES)) %>%
  column_to_rownames("bcr_patient_barcode")

#GLM: Geneset Genes
glm10 <- glm(response_binary~.,
    data = modeldata, family = "binomial")

  summary(glm10)$coef

#Reduce model stepwise by AIC
fit <- glm10 %>% step()

  summary(fit)$coef %>%  class()
  

# smry <- summary(fit)$coef
# fint <- confint(fit)
# OR <- exp(smry[,"Estimate"])
# lower_bound <- exp(fint[,1])
# upper_bound <- exp(fint[,2])

#Evaluate Model
#Area Under The Curve ####
pred <- predict(fit,type="response",data=modeldata)
roc(modeldata$response_binary~pred) %>% 
plot.roc(.,print.thres = T,print.auc = T)

save.image("capstone.RData")

####PLOTTING####
plot(PCA)
plot(norm_PCA)
min(pval$p, na.rm= TRUE) %>% print()

ggplot(pval, aes(x=p)) +
  geom_histogram(alpha=0.4) +
  labs(x=NULL, y=NULL, title="P-value")
