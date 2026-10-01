# Setup- install packages and format data from cvs file-------------------------

# Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install WGCNA's Bioconductor dependencies
BiocManager::install(c("impute", "GO.db", "preprocessCore"))

# Install and open WGCNA
install.packages("WGCNA")

library(WGCNA)

# Read in your raw expression file
microarraytable <- read.csv("GSE35182_ReadyToDeDupe.csv", header = TRUE)

View(microarraytable)

# Remove rows with ambiguous gene symbols (multiple genes joined by "///") 
microarraytable2 <- microarraytable[!grepl("///", microarraytable$Gene.Symbol), ]

# Build expression matrix: probe x samples
expr_matrix <- as.matrix(microarraytable2[, c("PBS_1","PBS_2","PBS_3","PBS_4","PBS_5","PBS_6",
                                    "CVB3_1","CVB3_2","CVB3_3","CVB3_4","CVB3_5","CVB3_6")])

# Row names = probe IDs 
rownames(expr_matrix) <- microarraytable2$ID_REF

# Gene symbol for each row
gene_symbols <- microarraytable2$Gene.Symbol

# Duplicate removal ------------------------------------------------------------

# For genes with multiple probes, keep the probe with the highest mean expression
result <- collapseRows(datET = expr_matrix,
                       rowGroup = gene_symbols,
                       rowID = microarraytable2$ID_REF,
                       method = "MaxMean")   

# The collapsed expression matrix (one row per gene now)
collapsed_expr <- result$datETcollapsed

# Which probe ID had the highest mean expression for each gene
selected_probes <- result$group2row

# For doccumentation purposes
write.csv(collapsed_expr, "collapsed_expression.csv", row.names = TRUE)

# Rebuild a clean output table: ID + Gene Symbol + expression values-----------

final_data <- data.frame(
  ID_REF = selected_probes[, "selectedRowID"],
  Gene.Symbol = rownames(collapsed_expr),
  collapsed_expr,
  row.names = NULL
)

# Reorder columns: ID, Gene Symbol, then samples
final_data <- final_data[, c("ID_REF", "Gene.Symbol", 
                             "PBS_1","PBS_2","PBS_3","PBS_4","PBS_5","PBS_6",
                             "CVB3_1","CVB3_2","CVB3_3","CVB3_4","CVB3_5","CVB3_6")]

write.csv(final_data, "GSE35182_DeDuped.csv", row.names = FALSE)

# Verification: confirm ID pairs to the correct gene symbol --------------

all(check$Gene.Symbol == check$Original.Symbol)
final_data$Original.Symbol <- microarraytable2$Gene.Symbol[
  match(final_data$ID_REF, microarraytable2$ID_REF)
]

all(final_data$Gene.Symbol == final_data$Original.Symbol)
