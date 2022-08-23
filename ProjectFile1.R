setwd("~/")
pacman::p_load(tidyverse, janitor)
rm(list=ls())



##read tgca files

tcga_case_data<-read_excel("TCGA drug response_Rboot.xlsx")
tcga_coad<-read.csv("TCGA_COAD_TMM_p400_common_Rboot.csv")
tcga_read<-read.csv("TCGA_READ_TMM_p400_common_Rboot.csv")
tcga_case_data$measure_of_response


##cataloge according to response
## cateogries of response : unique(tcga_case_data$measure_of_response)
diag<-unique(tcga_case_data$measure_of_response)

resp <- tcga_case_data %>% select(Cancer,bcr_patient_barcode, drug_name, measure_of_response)
resp %>%  
  filter(Cancer=="Colon adenocarcinoma (COAD)", drug_name=="Fluorouracil") %>%
  tabyl(measure_of_response)
  
resp_data<-resp %>%  
  filter(Cancer=="Colon adenocarcinoma (COAD)") %>%
  mutate(value=1) %>% 
  distinct() %>% 
  pivot_wider(names_from = drug_name, values_from = value)
  
FU_resp_data<-resp_data %>%  
  filter(Fluorouracil==1)

patient_id<-colnames(tcga_coad)[-1]
patient_id<-gsub("\\.", "-", patient_id)
patient_id<-str_split_fixed(patient_id, "-", 4)[,-4]
patient_id<-paste0(patient_id[,1],"-", patient_id[,2],"-", patient_id[,3] )

# work in progress:
# patient_id<-colnames(tcga_coad)[-1] %>% 
#   gsub("\\.", "-", .) %>% 
#   str_split_fixed(., "-", 4)[,-4]%>% 
#   paste0([,1],"-", [,2],"-", [,3] )

#patient ids common between response and expression dataset
coad_patient_id<-unique(intersect(patient_id,FU_resp_data$bcr_patient_barcode))


#remove patient without expression data
setdiff(FU_resp_data$bcr_patient_barcode,coad_patient_id)

FU_resp_data2<-FU_resp_data %>% filter(!bcr_patient_barcode%in%setdiff(FU_resp_data$bcr_patient_barcode,coad_patient_id))

FU_resp_data2 %>% tabyl(measure_of_response)

colnames(tcga_coad)<-c("gene_id",patient_id)

tcga_coad_FU<-tcga_coad %>% select(c("gene_id", coad_patient_id))

tcga_coad_FU_tp<-tcga_coad_FU %>% t()
