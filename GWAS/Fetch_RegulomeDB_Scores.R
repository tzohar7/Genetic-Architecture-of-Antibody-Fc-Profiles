########################################################
# GET REGULOMEDB SCORES
########################################################

library(org.Hs.eg.db)
library("data.table")
library("readr")

# CHANGE tbl file if needed
filename <- "~data/results_table.txt"
df <- read.table(filename, sep = '\t', header = T)

regulomefile <- "~data/databases/RegulomeDB.tsv"
reg <- fread(regulomefile, sep='\t')
reg <- as.data.frame(reg)

df[names(reg)[5:15]] <- ""

idx <- which(reg$rsid %in% df$rsID)
reg <- reg[idx,]

#df[i,c("Ref", "Alt")]

for (i in 1:dim(df)[1]) {
	regi <- reg[reg$rsid == df$rsID[i],]
	tru <- (regi[,c("ref", "alt")] == df[i,"Ref"]) + (regi[,c("ref", "alt")] == df[i,"Alt"])
	regi <- regi[apply(tru,1,sum) == 2,]
	if (dim(regi)[1] != 0) {
		df[i, names(reg)[5:15]] <- regi[,5:15]
}
}

# CHANGE tbl file if needed
write.table(df, "~data/results_table.txt", append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
