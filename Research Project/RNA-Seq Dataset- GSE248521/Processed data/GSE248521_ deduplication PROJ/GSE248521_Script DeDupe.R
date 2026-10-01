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
rnaseqtable <- read.csv("GSE248521_ReadyToDeDupe.csv", header = TRUE)

View(rnaseqtable)

# Remove rows with ambiguous gene symbols  
rnaseqtable2 <- rnaseqtable[!grepl("///|#or#", rnaseqtable$gene_symbol), ]

# Build expression matrix: gene x samples
expr_matrix <- as.matrix(rnaseqtable2[, c("sham_1","sham_2","sham_3","sham_4","sham_5","sham_6",
                                          "CVB3_1","CVB3_2","CVB3_3","CVB3_4","CVB3_5","CVB3_6")])

# Row names = original IDs 
rownames(expr_matrix) <- rnaseqtable2$id

# Gene symbol for each row
gene_symbols <- rnaseqtable2$gene_symbol

# Duplicate removal ------------------------------------------------------------

# For genes with multiple IDs, keep the one with the highest mean expression
result <- collapseRows(datET = expr_matrix,
                       rowGroup = gene_symbols,
                       rowID = rnaseqtable2$id,
                       method = "MaxMean")   

# The collapsed expression matrix (one row per gene now)
collapsed_expr <- result$datETcollapsed

# Which ID had the highest mean expression for each gene
selected_ids <- result$group2row

# For doccumentation purposes
write.csv(collapsed_expr, "collapsed_expression_rnaseq.csv", row.names = TRUE)

# Rebuild a clean output table: ID + Gene Symbol + expression values-----------

final_data <- data.frame(
  id = selected_ids[, "selectedRowID"],
  gene_symbol = rownames(collapsed_expr),
  collapsed_expr,
  row.names = NULL
)

# Reorder columns: ID, Gene Symbol, then samples
final_data <- final_data[, c("id", "gene_symbol", 
                             "sham_1","sham_2","sham_3","sham_4","sham_5","sham_6",
                             "CVB3_1","CVB3_2","CVB3_3","CVB3_4","CVB3_5","CVB3_6")]

write.csv(final_data, "GSE248521_DeDuped.csv", row.names = FALSE)

# Verification: confirm ID pairs to the correct gene symbol --------------

final_data$Original.Symbol <- rnaseqtable2$gene_symbol[
  match(final_data$id, rnaseqtable2$id)
]

all(final_data$gene_symbol == final_data$Original.Symbol)
