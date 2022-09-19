########################################################
# MAKE TO AVINPUT FILE IN R 
########################################################

filename <- "~data/results/results_table.txt"
df_avout <- read.table(filename, sep = "\t", header = TRUE)

df_avout <- df_avout[,c("CHR","BP","BP","Ref","Alt","rsID")]
#df_avout <- df_avout[,c("CHR","BP","BP","A1","A2","rsID","P", "BETA", "SE", "Phenotype", "Genes")]

write.table(df_avout, "~data/results/results_table.avinput", append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
