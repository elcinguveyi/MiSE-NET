#for GWAS
gwas_catalog = read.table(file = 'forGWAS/GWAS_Catalog/gwas-association-downloaded_2024-09-18-EFO_0001645-withChildTraits.tsv', sep = '\t', header = TRUE, fill=TRUE)
gwas_cardiogram = read.table(file = '0_previous/fromCardiogram/cardiogram_gwas_results/CARDIoGRAM_GWAS_RESULTS.txt', sep = '\t', header = TRUE)

#for eQTL
