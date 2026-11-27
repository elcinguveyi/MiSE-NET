#for GWAS
gwas_catalog = read.table(file = 'forGWAS/GWAS_Catalog/gwas-association-downloaded_2024-09-18-EFO_0001645-withChildTraits.tsv', sep = '\t', header = TRUE, fill=TRUE)
gwas_cardiogram = read.table(file = 'forGWAS/fromCardiogram/cardiogram_gwas_results/CARDIoGRAM_GWAS_RESULTS.txt', sep = '\t', header = TRUE)

#for eQTL
#GTEx
tissues = list("Artery_Aorta","Artery_Coronary","Artery_Tibial","Heart_Atrial_Appendage","Heart_Left_Ventricle") #for heart tissue
tissues = list("Whole_Blood") #for blood tissue


for (item in tissues){
  R.utils::gunzip(glue::glue('forEQTL/GTEX/GTEx_Analysis_v8_eQTL/{item}.v8.egenes.txt.gz'))  
  eqtl = read.table(file = glue::glue('forEQTL/GTEX/GTEx_Analysis_v8_eQTL/{item}.v8.egenes.txt'), sep = '\t', header = TRUE)
  name <- paste("eqtl", item, sep = "_")
  assign(name,  as.data.frame(eqtl))
}

library(tidyverse)
eqtl_Artery_Aorta = eqtl_Artery_Aorta %>% add_column(tissue_name = "Artery_Aorta")
eqtl_Artery_Coronary = eqtl_Artery_Coronary %>% add_column(tissue_name = "Artery_Coronary")
eqtl_Artery_Tibial = eqtl_Artery_Tibial %>% add_column(tissue_name = "Artery_Tibial")
eqtl_Heart_Atrial_Appendage = eqtl_Heart_Atrial_Appendage %>% add_column(tissue_name = "Heart_Atrial_Appendage")
eqtl_Heart_Left_Ventricle = eqtl_Heart_Left_Ventricle %>% add_column(tissue_name = "Heart_Left_Ventricle")
eqtl_heart_gtex = rbind(eqtl_Artery_Aorta,eqtl_Artery_Coronary,eqtl_Artery_Tibial,eqtl_Heart_Atrial_Appendage,eqtl_Heart_Left_Ventricle) 
write.csv(eqtl_heart_gtex, "eqtl_heart_gtex.csv") 
eqtl_heart_gtex = read.table("eqtl_heart_gtex.csv", sep = ',', header = TRUE, row.names = 1)

write.csv(eqtl_Whole_Blood, "eqtl_blood_gtex.csv") 
eqtl_blood_gtex = read.table("eqtl_blood_gtex.csv", sep = ',', header = TRUE, row.names = 1)

#STARNET
#all modules
library(dplyr)
folder_path <- "~/Desktop/GWAS/forEQTL/STARNET/all/"
file_list <- list.files(path = folder_path, pattern = "*.csv", full.names = TRUE)
df_list <- lapply(file_list, function(file) {
  df <- read.csv(file, stringsAsFactors = FALSE)  # CSV'yi oku
  df <- df %>% filter(tissue == "AOR")  
  return(df)
})
final_df <- bind_rows(df_list)

final_df <- final_df %>% filter(!is.na(hgnc_symbol))

duplicate_rows <- final_df %>% 
  group_by(SNP, hgnc_symbol) %>% 
  filter(n() > 1)
print(duplicate_rows) # not found

eqtl_starnet <- final_df
write.csv(eqtl_starnet, "eqtl_*tissueName_starnet.csv") #tissueName: heart or blood
eqtl_heart_starnet = read.table("eqtl_heart_starnet.csv", sep = ',', header = TRUE, row.names = 1)
eqtl_blood_starnet = read.table("eqtl_blood_starnet.csv", sep = ',', header = TRUE, row.names = 1)


#FILTER GWAS
pval = 0.00001
filtered_gwas_cardiogram = subset(gwas_cardiogram, pvalue < pval) 
filtered_gwas_cardiogram = filtered_gwas_cardiogram[, c("SNP", "pvalue")]
colnames(filtered_gwas_cardiogram) = c("SNP", "PVALUE")
filtered_gwas_catalog = subset(gwas_catalog, P.VALUE < pval) 
filtered_gwas_catalog = filtered_gwas_catalog[, c("SNPS","P.VALUE")]
colnames(filtered_gwas_catalog) = c("SNP", "PVALUE")

all_gwas = rbind(filtered_gwas_cardiogram,filtered_gwas_catalog) 
length(unique(all_gwas$SNP)) 
rownames(all_gwas) = NULL
all_gwas <- all_gwas[order(all_gwas$SNP), ] #first rows ch:.. chromosome locations
all_gwas <- subset(all_gwas, grepl("rs", SNP)) #filter including 'rs' 
write.csv(all_gwas, "all_gwas_10-5.csv")
unique_SNP_list = unique(all_gwas$SNP)
unique_SNP_list = as.data.frame(unique_SNP_list)
write.csv(unique_SNP_list, "unique_SNP_list_10-5.csv")


all_gwas = rbind(filtered_gwas_cardiogram,filtered_gwas_catalog) 
count(unique(all_gwas["SNP"])) 

tekrarliSNP <- all_gwas[all_gwas$SNP %in% all_gwas$SNP[duplicated(all_gwas$SNP)], ]
result <- anti_join(all_gwas, tekrarliSNP) 
count(unique(all_gwas["SNP"])) - count(unique(tekrarliSNP["SNP"]))
tekrarliSNP <- tekrarliSNP %>% group_by(SNP) %>%  filter(PVALUE == min(PVALUE)) %>% slice(1) %>% ungroup() 
count(unique(tekrarliSNP["SNP"])) 
all_gwas2 = rbind(result, tekrarliSNP)
count(unique(all_gwas2["SNP"])) # all unique


#FILTER EQTL
#selection
eqtl_starnet = eqtl_heart_starnet 
eqtl_gtex = eqtl_heart_gtex 

filtered_eqtl_starnet = eqtl_starnet[, c("SNP","convert","p.value")] 
colnames(filtered_eqtl_starnet) = c("SNP", "GENE", "PVALUE")
filtered_eqtl_starnet <- filtered_eqtl_starnet %>% filter(!is.na(GENE) & str_trim(GENE) != "") 
filtered_eqtl_gtex = eqtl_gtex[,c("rs_id_dbSNP151_GRCh38p7","gene_name","pval_nominal")]
colnames(filtered_eqtl_gtex) = c("SNP", "GENE", "PVALUE")

all_eqtl = rbind(filtered_eqtl_gtex, filtered_eqtl_starnet)
length(unique(all_eqtl$GENE)) 
length(unique(all_eqtl$SNP)) 
write.csv(all_eqtl, "all_eqtl_liver.csv")

all_gwas = read.table("all_gwas_10-5.csv", sep = ',', header = TRUE, row.names = 1) 
all_eqtl = read.table("all_eqtl_blood.csv", sep = ',', header = TRUE, row.names = 1)

#for GEO 

