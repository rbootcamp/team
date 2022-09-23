setwd("~/capstone")
pacman::p_load(tidyverse, janitor, readxl, ggfortify, ggrepel)
rm(list=setdiff(ls(),'v'))
theme_set( theme_bw() )


#Read TGCA files ####
tcga_case_data<-read_csv("TCGA drug response_Rboot_2.csv")
tcga_coad<-read_csv("TCGA_COAD_TMM_p400_common_Rboot.csv")
tcga_read<-read_csv("TCGA_READ_TMM_p400_common_Rboot.csv")
tcga_coad_read<-merge(tcga_coad,tcga_read,by="gene_id")
rm(tcga_coad, tcga_read)

#Catalogue according to response ####
#cateogries of response : unique(tcga_case_data$measure_of_response)
unique(tcga_case_data$measure_of_response)

tcga_case_data %>%
  select(Cancer,bcr_patient_barcode, drug_name, measure_of_response) %>%  
  filter(drug_name=="Fluorouracil") %>%
  tabyl(measure_of_response)

FU_resp_data <- tcga_case_data %>% select(Cancer,bcr_patient_barcode, drug_name, measure_of_response) %>%  
  mutate(value=1) %>% 
  distinct() %>% 
  pivot_wider(names_from = drug_name, values_from = value) %>%  
  filter(Fluorouracil==1)

patient_id<-colnames(tcga_coad_read)[-1] %>%
  gsub("\\.", "-", .)

patient_id<-str_split_fixed(patient_id, "-", 4)[,-4]
patient_id<-paste0(patient_id[,1],"-", patient_id[,2],"-", patient_id[,3] )

# work in progress:
# patient_id<-colnames(tcga_coad)[-1] %>% 
#   gsub("\\.", "-", .) %>% 
#   str_split_fixed(., "-", 4)[,-4]%>% 
#   paste0([,1],"-", [,2],"-", [,3] )

#Patient ids common between response and expression dataset ####
coad_read_patient_id<-unique(intersect(patient_id,
                                       FU_resp_data$bcr_patient_barcode))

#Remove patient without expression data ####
setdiff(FU_resp_data$bcr_patient_barcode,
        coad_read_patient_id)

FU_resp_data2<-FU_resp_data %>%
  filter(!bcr_patient_barcode%in%
           setdiff(FU_resp_data$bcr_patient_barcode,
                   coad_read_patient_id))

rm(tcga_case_data, FU_resp_data)

FU_resp_data2 %>% tabyl(measure_of_response)

colnames(tcga_coad_read)<-c("gene_id",patient_id)

tcga_coad_read_FU_tp<-tcga_coad_read %>%
  select(c("gene_id", patient_id)) %>%
  column_to_rownames("gene_id") %>%
  t() %>%
  data.frame() %>%
  rownames_to_column("bcr_patient_barcode")

rm(patient_id, tcga_coad_read, coad_read_patient_id)

#Merging gene expression and drug response table ####

FU_resp_data2 %>% dim()
tcga_coad_read_FU_tp %>% dim()

intersect(FU_resp_data2$bcr_patient_barcode,
          tcga_coad_read_FU_tp$bcr_patient_barcode)

# load("comb.rda")

comb <- inner_join(FU_resp_data2, tcga_coad_read_FU_tp) %>%
  mutate(response_status = ifelse(measure_of_response %in% c("Partial Response", "Complete Response"),
                                  "Responder", "Nonresponder"), .after = measure_of_response) %>%
  mutate(response_binary = recode(response_status,
                                  "Responder" = 1, "Nonresponder" = 0), .after = response_status)
comb %>% dim()

save(comb, file="comb.rda")
rm(FU_resp_data2, tcga_coad_read_FU_tp)

#make vector of gene IDs
gid <- comb %>%
  select(starts_with("ENSG")) %>%
  colnames()

#NORMALISATION FUNCTION#
cal_log2fc <- function(df) {
  lg <- log2(df+1)
  med <- log2(median(df+1))
  output <- lg-med
  return(output)
}
#FOR LOOP FOR NORMALISATION ####
# load("nrm.rda")

# ap <- sapply(comb[,gid], function(x) cal_log2fc(x) )
# ap <- sapply(comb[,gid], cal_log2fc )

nrm <- comb

#sapply method (faster)
nrm[,gid] <- sapply(comb[,gid], cal_log2fc )
save(nrm, file="nrm.rda")

#for loop method...
for(i in gid){
  nrm[,i] <- comb %>%
    pull(i) %>%
    cal_log2fc()
}
save(nrm, file="nrm.rda")

#FOR LOOP FOR WILCOX-TEST ####
# load("pval.rda")

#make dataframe for storing p-values
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

#Data exploration and Visualisation ####

min(pval$p, na.rm= TRUE)

pval %>% filter(p < 0.05/54098)

ggplot(pval, aes(x=p)) +
  geom_histogram(alpha=0.4) +
  labs(x=NULL, y=NULL, title="P-value")


# PRINCIPAL COMPONENT ANALYSIS #

# calculate variance of each column in expression data
# load("variance_table.rda")
v <- comb %>%
  select(all_of(gid)) %>% 
  var(use="pairwise.complete.obs") %>% 
  diag()

save(v, file="variance_table.rda")

#list of genes with non-zero variance
keep <- which(v!=0) %>%
  names()

#PCA (RAW DATA)
pc <- prcomp(comb %>% select(all_of(keep)),
             scale = TRUE)

summary(pc)

autoplot(pc, data=comb, col="response_status")

#PCA (NORMALISED DATA)
pc_nrm <- prcomp(nrm %>% select(all_of(keep)),
                 scale = TRUE)

summary(pc_nrm)

autoplot(pc_nrm, data=nrm, col="response_status")

# Test GLM #

head(comb$response_binary)
fit <- glm(response_binary~as.numeric(ENSG00000119396)+as.numeric(ENSG00000119396),data = comb, family = "binomial")
summary(fit)
