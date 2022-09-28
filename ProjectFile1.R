setwd("~/capstone")
pacman::p_load(tidyverse, janitor, readxl, ggfortify, ggrepel, DataExplorer, ClassDiscovery, pROC)
rm(list=ls())
theme_set( theme_bw() )

#Read TGCA files ####
tcga_case_data<-read_csv("TCGA drug response_Rboot_2.csv")
tcga_coad<-read_csv("TCGA_COAD_TMM_p400_common_Rboot.csv")
tcga_read<-read_csv("TCGA_READ_TMM_p400_common_Rboot.csv")
tcga_coad_read<-merge(tcga_coad,tcga_read,by="gene_id")
rm(tcga_coad, tcga_read)

#Catalogue according to response ####

#Print response category names
unique(tcga_case_data$measure_of_response) %>%
  print()

#Print a table of responses (Fluorouracil recipients)
tcga_case_data %>%
  select(Cancer,bcr_patient_barcode, drug_name, measure_of_response) %>%  
  filter(drug_name=="Fluorouracil") %>%
  tabyl(measure_of_response)

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

#Remove patients without expression data ####

#Make vector of patient ids common between response and expression datasets
common_id <- intersect(patient_id,
                       fcil_case_data$bcr_patient_barcode) %>%
  unique()

#Identify patients without expression data
setdiff(fcil_case_data$bcr_patient_barcode,
        common_id) %>%
  print()

#Remove patients without expression data
response_data <- fcil_case_data %>%
  filter(!bcr_patient_barcode%in%
           setdiff(fcil_case_data$bcr_patient_barcode,
                   common_id))

#Response table - Fluorouracil patients with expression data#
response_data %>%
  tabyl(measure_of_response)

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

#Merge gene expression and drug response table ####
# load("comb.rda")

#Print ncol and nrow of each dataframe
response_data %>% dim()
gene_matrix %>% dim()

#Construct combined dataset
comb <- inner_join(response_data, gene_matrix) %>%
  select( !which(colSums(.[,sapply(.,is.numeric)])==0) %>% names() ) %>% #filter columns with no reads
  mutate(response_status = ifelse(measure_of_response %in% c("Partial Response", "Complete Response"),
                                  "Responder", "Nonresponder"), .after = measure_of_response) %>%
  mutate(response_binary = recode(response_status, #add columns for binomial response status
                                  "Responder" = 1, "Nonresponder" = 0), .after = response_status)
#New dimensions
comb %>% dim()

save(comb, file="comb.rda")

rm(response_data, gene_matrix)

#Make vector of gene IDs
gid <- comb %>%
  select(starts_with("ENSG")) %>%
  colnames()

#Normalise####
#NORMALISATION FUNCTION#
cal_log2fc <- function(df) {
  lg <- log2(df+1)
  med <- log2(median(df+1))
  output <- lg-med
  return(output)
}
#sapply for normalisation####
# load("nrm.rda")

nrm <- comb

nrm[,gid] <- sapply(comb[,gid], cal_log2fc )
save(nrm, file="nrm.rda")

#Data exploration and Visualisation ####

# PRINCIPAL COMPONENT ANALYSIS ####

#PCA (RAW DATA)
pc <- prcomp(comb %>% select(all_of(gid)),
             scale = TRUE)

summary(pc)

autoplot(pc, data=comb, col="response_status")

#PCA (NORMALISED DATA)
pc_nrm <- prcomp(nrm %>% select(all_of(gid)),
                 scale = TRUE)

summary(pc_nrm)

autoplot(pc_nrm, data=nrm, col="response_status")

#FOR LOOP FOR WILCOX-TEST ####
# load("pval.rda")

#Make dataframe for storing p-values
pval <- gid %>%
  data.frame(gene=.) %>% 
  add_column(p=NA) %>% 
  column_to_rownames('gene')
head(pval)

# Wilcox-test
for(i in gid){
  pval[i,"p"] <- wilcox.test(as.numeric(as.matrix(nrm[,i])) ~
                               nrm$response_status,
                               na.action = na.omit)$p.value
}

save(pval, file="pval.rda")
View(pval)

#Explore significant genes ####

ggplot(pval, aes(x=p)) +
  geom_histogram(alpha=0.4) +
  labs(x=NULL, y=NULL, title="P-value")

min(pval$p, na.rm= TRUE) %>% print()

pval %>% filter(p < 0.05/nrow(pval)) %>% print()

print(pval_top11 <- pval %>% filter(p < 1e-5))

#Make vector of top 11 genes
top11 <- pval_top11 %>% rownames()

#HEATMAP####

#Subset top 11 significant genes
mapdata <- nrm %>%
  select(bcr_patient_barcode, all_of(top11)) %>%
  column_to_rownames("bcr_patient_barcode")

#Correlation Matrix
plot_correlation(mapdata %>%
                   as.matrix())

#The most highly correlated genes
cor.test(mapdata$ENSG00000250956, mapdata$ENSG00000230355)

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

#Model top11 genes
nrm %>%
  select(bcr_patient_barcode, response_binary, all_of(top11)) %>%
  column_to_rownames("bcr_patient_barcode") %>%
  glm(response_binary~.,
      data = ., family = "binomial") %>%
  summary()

#Gene "ENSG00000207395" has a very high Std. Error (622.59703) -- removed from further analysis
modeldata <- nrm %>%
  select(bcr_patient_barcode, response_binary, all_of(top11[top11!="ENSG00000207395"])) %>%
  column_to_rownames("bcr_patient_barcode")

#GLM: 10 Genes
glm10 <- glm(response_binary~.,
    data = modeldata, family = "binomial")

  summary(glm10)

#Reduce model stepwise by AIC
fit <- glm10 %>% step()

  summary(fit)$coef

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
