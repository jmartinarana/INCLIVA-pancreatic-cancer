library("dndscv")

#data("dataset_simbreast", package="dndscv")

########################### tumor ######################

tumor_input = "/home/mibarrola/UMP444_exomas_pancreas_1108/dnds/input_dnds_tejido.txt"

GRCh38_rda <- "/home/szuniga/UB17_exomas_CabezaCuello/scripts/dnds/GRCh38_refcds.rda"
out_folder <- "/home/mibarrola/UMP444_exomas_pancreas_1108/dnds/results/tejido/"
#ci = geneci(dndsout2)
tumor = read.csv(tumor_input, header = TRUE, sep = "\t", encoding = "ASCII", stringsAsFactors = FALSE)

muestras = unique(tumor$sampleID)
for ( mu in muestras) {
  p8 = tumor[which(tumor$sampleID == mu),]
  print (mu)
  dndsout2 = dndscv(p8, refdb=GRCh38_rda,
                    max_muts_per_gene_per_sample = 1000, max_coding_muts_per_sample = 10000, outmats = T, cv = NULL)
  sel_cv = dndsout2$sel_cv
  signif_genes = sel_cv[sel_cv$qglobal_cv<0.05, c("gene_name","wmis_cv","wnon_cv","wspl_cv","wind_cv", "qglobal_cv")]
  write.table(signif_genes, paste0(out_folder, mu, "signif_genes_tumor.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)
  global = dndsout2$globaldnds
  write.table(global, paste0(out_folder, mu, "global_tumor.txt"), sep = "\t", col.names = TRUE, row.names = TRUE)
}


################# plasmas ###########

plasma_input = "/home/szuniga/UB17_exomas_CabezaCuello/results/otros_analisis/dnds/input_plasma_dnds.txt"
tumor = read.csv(plasma_input, header = TRUE, sep = "\t", encoding = "ASCII", stringsAsFactors = FALSE)


GRCh38_rda <- "/home/szuniga/UB17_exomas_CabezaCuello/scripts/dnds/GRCh38_refcds.rda"
out_folder <- "/home/mibarrola/UMP444_exomas_pancreas_1108/dnds/results/plasma/"
#ci = geneci(dndsout2)

muestras = unique(tumor$sampleID)
for ( mu in muestras) {
  p8 = tumor[which(tumor$sampleID == mu),]
  print (mu)
  dndsout2 = dndscv(p8, refdb=GRCh38_rda,
                    max_muts_per_gene_per_sample = 1000, max_coding_muts_per_sample = 10000, outmats = T, cv = NULL)
  sel_cv = dndsout2$sel_cv
  signif_genes = sel_cv[sel_cv$qglobal_cv<0.05, c("gene_name","wmis_cv","wnon_cv","wspl_cv","wind_cv", "qglobal_cv")]
  write.table(signif_genes, paste0(out_folder, mu, "signif_genes_tumor.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)
  global = dndsout2$globaldnds
  write.table(global, paste0(out_folder, mu, "global_tumor.txt"), sep = "\t", col.names = TRUE, row.names = TRUE)
}



######## plasmas1 ###########

tumor = read.csv("/home/szuniga/UB17_exomas_CabezaCuello/results/otros_analisis/dnds/input_plasma1_dnds.vcf", header = TRUE, sep = "\t", encoding = "ASCII", stringsAsFactors = FALSE)


GRCh38_rda <- "/home/szuniga/UB17_exomas_CabezaCuello/scripts/dnds/GRCh38_refcds.rda"
out_folder <- "/home/szuniga/UB17_exomas_CabezaCuello/results/otros_analisis/dnds/results/plasma1/"
#ci = geneci(dndsout2)

muestras = unique(tumor$sampleID)
for ( mu in muestras) {
  p8 = tumor[which(tumor$sampleID == mu),]
  print (mu)
  dndsout2 = dndscv(p8, refdb=GRCh38_rda,
                    max_muts_per_gene_per_sample = 1000, max_coding_muts_per_sample = 10000, outmats = T, cv = NULL)
  sel_cv = dndsout2$sel_cv
  signif_genes = sel_cv[sel_cv$qglobal_cv<0.05, c("gene_name","wmis_cv","wnon_cv","wspl_cv","wind_cv", "qglobal_cv")]
  write.table(signif_genes, paste0(out_folder, mu, "signif_genes_tumor.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)
  global = dndsout2$globaldnds
  write.table(global, paste0(out_folder, mu, "global_tumor.txt"), sep = "\t", col.names = TRUE, row.names = TRUE)
}
