setwd("~/capstone")
pacman::p_load(tidyverse, janitor, readxl, ggplot2)
rm(list=ls())
theme_set( theme_bw() )


##read tgca files

tcga_case_data<-read_csv("TCGA drug response_Rboot_2.csv")
tcga_coad<-read_csv("TCGA_COAD_TMM_p400_common_Rboot.csv")
tcga_read<-read_csv("TCGA_READ_TMM_p400_common_Rboot.csv")
tcga_coad_read<-merge(tcga_coad,tcga_read,by="gene_id")
rm(tcga_coad, tcga_read)


##cataloge according to response
## cateogries of response : unique(tcga_case_data$measure_of_response)
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

#patient ids common between response and expression dataset
coad_read_patient_id<-unique(intersect(patient_id,
                                       FU_resp_data$bcr_patient_barcode))


#remove patient without expression data
setdiff(FU_resp_data$bcr_patient_barcode,
        coad_read_patient_id)

FU_resp_data2<-FU_resp_data %>%
  filter(!bcr_patient_barcode%in%
           setdiff(FU_resp_data$bcr_patient_barcode,
                   coad_read_patient_id))

rm(FU_resp_data)

FU_resp_data2 %>% tabyl(measure_of_response)

colnames(tcga_coad_read)<-c("gene_id",patient_id)

tcga_coad_read_FU_tp<-tcga_coad_read %>%
  select(c("gene_id", patient_id)) %>%
  column_to_rownames("gene_id") %>%
  t() %>%
  data.frame() %>%
  rownames_to_column("bcr_patient_barcode")

rm(patient_id, tcga_coad_read)


#Merging gene expression and drug response table


FU_resp_data2 %>% dim()
tcga_coad_read_FU_tp %>% dim()

intersect(FU_resp_data2$bcr_patient_barcode,
          tcga_coad_read_FU_tp$bcr_patient_barcode)

comb <- inner_join(FU_resp_data2, tcga_coad_read_FU_tp) %>%
  mutate(response_status = recode(measure_of_response, "Clinical Progressive Disease"= "Nonresponder", 
                            "Stable Disease" = "Nonresponder", "Complete Response" = "Responder", 
                            "Partial Response" = "Responder"), .after = measure_of_response) %>%
  mutate(response_binary = ifelse(response_status %in% "Responder", 1, 0 ), .after=response_status)


comb %>% dim()

save(comb, file="comb.rda")
# load("comb.rda")

#Statistical test to identify genes that are differential between responder and non-responder?


# colnames(d[10:20])
# 
# 
# is.na(d)
# 
# getOption("na.action")
# 
# wilcox.test(as.numeric(ENSG00000000419) ~ response2, data=d, na.action = na.omit) %>% names()
# 
# 
# ENSG00000000419.p <- wilcox.test(as.numeric(ENSG00000000419) ~ response2, data=d, na.action = na.omit)$p.value
# 


#make dataframe for storing p-values

gid <- comb %>%
  select(starts_with("ENSG")) %>%
  colnames() %>%
  data.frame(gene=.) %>% 
  add_column(p=NA) %>% 
  column_to_rownames('gene')
head(gid)


# FOR LOOP FOR WILCOX-TEST #
  
  
for(i in rownames(gid)){
  gid[i,"p"] <- wilcox.test(as.numeric(as.matrix(comb[,i])) ~ comb$response_status, na.action = na.omit)$p.value
}

save(gid, file="gid.rda")
#load("gid.rda")

# Data exploration #

min(gid$p, na.rm= TRUE)

gid %>% filter(p < 0.05/54098)

ggplot(gid, aes(x=p)) +
  geom_histogram(alpha=0.4) +
  labs(x=NULL, y=NULL, title="P-value")

# Test GLM #

head(comb$response_binary)
fit <- glm(response_binary~as.numeric(ENSG00000119396)+as.numeric(ENSG00000119396),data = comb, family = "binomial")
summary(fit)
