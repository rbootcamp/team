setwd("~/capstone")
pacman::p_load(tidyverse, janitor, ggfortify, DataExplorer, pheatmap, pROC, Hmisc)
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
#load("comb.rda")

#Print ncol and nrow of each dataframe
response_data %>% dim()
gene_matrix %>% dim()

#Construct combined dataset
comb <- inner_join(response_data, gene_matrix) %>%
  select( !which(colSums(.[,sapply(.,is.numeric)])==0) %>% names() ) %>% #filter columns with no reads
  mutate(response_status = ifelse(measure_of_response %in% c("Partial Response", "Complete Response"),
                                  "Responder", "Nonresponder"), .after = measure_of_response) %>%
  mutate(response_binary = recode(response_status, #add columns for binomial response status
                                  "Responder" = 1, "Nonresponder" = 0), .after = response_status) %>%
  column_to_rownames("bcr_patient_barcode")

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
  geom_histogram(alpha=0.4, color="black", fill="lightblue") +
  labs(x="P-value", y="Frequency", 
       title="Frequency of genes with various P-value associated with 5-FU response")


min(pval$p, na.rm= TRUE) %>% print()

pval %>% filter(p < 0.05/nrow(pval)) %>% print()

print(pval_top11 <- pval %>% filter(p < 1e-5))


#Make vector of top 11 genes
top11 <- pval_top11 %>% rownames()

#Benjamin Hochberg correction
BHpval<-data.frame(rownames(pval)
                   [which((p.adjust(pval$p, method = "BH"))<=0.05)],
                   p.adjust(pval$p, method = "BH")
                   [which((p.adjust(pval$p, method = "BH"))<=0.05)])

#Make vector of BH genes
BHgenes<-BHpval[,1]
BHgenes<-BHgenes[1:16]

#HEATMAP####

#Subset top 11 significant genes
mapdata <- nrm %>%
  select(all_of(top11))

#Correlation Matrix
plot_correlation(mapdata %>%
                   as.matrix())

cor_matrix <- rcorr(mapdata %>% as.matrix())$P #Matrix of P-values for each correlation
cor_matrix[lower.tri(cor_matrix,diag=TRUE)]=NA #Remove redundant values

print(cor_table <- #Print as table ordered by significance
        cor_matrix %>%
        as.table() %>%
        as.data.frame() %>%
        na.omit() %>%
        rename(p=Freq, gene1=Var1, gene2=Var2) %>% 
        arrange(p))

#The most highly correlated genes
cor.test(mapdata$ENSG00000250956, mapdata$ENSG00000230355)

#Table to annotate by response
response_table <- nrm %>%
  select(response_status)


#Heatmap
pheatmap(mapdata,
         annotation_row = response_table,
         show_rownames=F,
         fontsize=7)

# GLM ####

#Model BH genes
glm16 <- nrm %>%
  select(response_binary, all_of(BHgenes)) %>%
  glm(response_binary~.,
      data = ., family = "binomial")
  summary(glm16)

#List of 'good' genes without large Std. Error
smry <- summary(glm16)$coef %>%
  as.data.frame() %>%
  rename(std_error="Std. Error") %>%
  select(std_error) %>%
  filter(std_error<1e3) %>%
  rownames()
top9 <- smry[-1] #Remove "(Intercept)"
rm(glm16, smry)

#Genes with high Std. Error removed from further analysis
modeldata <- nrm %>%
  select(response_binary, all_of(top9))


#GLM: 9 Genes
glm9 <- glm(response_binary~.,
             data = modeldata, family = "binomial")

  summary(glm9)


#Reduce model stepwise by AIC
fit <- glm9 %>% step()

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
