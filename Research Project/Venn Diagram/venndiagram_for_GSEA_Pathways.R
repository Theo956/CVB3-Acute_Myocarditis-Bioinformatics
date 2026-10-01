# 3-way Venn diagram with pathway names, using base R 
# Note to reader: Run the WHOLE script at once (Ctrl+A, then click "Source")


# Data Import -----------------------------------------------------------------------------
file1 <- "GSE35182_Sig_Pathways.csv"
file2 <- "GSE248521_Sig_Pathways.csv"
file3 <- "GSE19496_Sig_Pathways.csv"

raw1 <- read.csv(file1, stringsAsFactors = FALSE)
raw2 <- read.csv(file2, stringsAsFactors = FALSE)
raw3 <- read.csv(file3, stringsAsFactors = FALSE)

# Disease-label pathways to exclude -----------------------------------------------------
disease_labels <- c(
  "Staphylococcus aureus infection", "Pertussis", "Leishmaniasis",
  "Primary immunodeficiency", "Systemic lupus erythematosus",
  "Graft-versus-host disease", "Inflammatory bowel disease", "Malaria",
  "Allograft rejection", "Legionellosis",
  "Intestinal immune network for IgA production", "Asthma",
  "Herpes simplex virus 1 infection", "Coronavirus disease-COVID-19",
  "Influenza A", "Measles", "Rheumatoid arthritis",
  "Autoimmune thyroid disease", "Type I diabetes mellitus",
  "Primary bile acid biosynthesis"
)

# Disease-label pathway removal ------------------------------------------------------------
get_filtered_pathways <- function(df, disease_labels) {
  col_idx <- grep("Pathways$", colnames(df))
  if (length(col_idx) == 0) {
    stop("Could not find a pathway-name column, check colnames(df)")
  }
  p <- trimws(df[[col_idx[1]]])
  p <- p[p != ""]
  p <- unique(p)
  p[!p %in% disease_labels]
}

A <- get_filtered_pathways(raw1, disease_labels)   
B <- get_filtered_pathways(raw2, disease_labels)   
C <- get_filtered_pathways(raw3, disease_labels)   

cat("Kept pathways per dataset:", length(A), length(B), length(C), "\n")

# Venn Diagram set up -----------------------------------------------------------------
reg <- list(
  A   = sort(setdiff(A, union(B, C))),
  B   = sort(setdiff(B, union(A, C))),
  C   = sort(setdiff(C, union(A, B))),
  AB  = sort(setdiff(intersect(A, B), C)),
  AC  = sort(setdiff(intersect(A, C), B)),
  BC  = sort(setdiff(intersect(B, C), A)),
  ABC = sort(Reduce(intersect, list(A, B, C)))
)
View(reg)

# Break up long pathway names onto separate lines with bullet points
make_label <- function(items, width) {
  if (length(items) == 0) return("(0)")
  wrapped <- sapply(items, function(x) paste(strwrap(x, width, initial = "\u2022 ", prefix = "   "), collapse = "\n"))
  paste0("(", length(items), ")\n", paste(wrapped, collapse = "\n"))
}

# Venn Diagram creation and save ------------------------------------------------------------------
r <- 3.4
cx <- c(A = -1.6, B = 1.6, C = 0)
cy <- c(A = 1.0,  B = 1.0, C = -1.7)
cols <- c(A = "#4C72B0", B = "#DD8452", C = "#55A868")

png("pathway_venn_R.png", width = 3000, height = 2400, res = 200)

par(mar = c(1, 1, 4, 1))
plot(NA, xlim = c(-7, 7), ylim = c(-6.3, 5.6), asp = 1,
     axes = FALSE, xlab = "", ylab = "",
     main = "Overlap of significant KEGG pathways across 3 datasets\n(disease-label duplicates excluded)")

for (k in c("A", "B", "C")) {
  symbols(cx[k], cy[k], circles = r, inches = FALSE, add = TRUE,
          bg = adjustcolor(cols[k], alpha.f = 0.35), fg = "black", lwd = 2)
}

# region text: x, y, key, wrap width
text(-3.4,  1.9, make_label(reg$A,   22), cex = 0.85)
text( 3.4,  1.9, make_label(reg$B,   22), cex = 0.85)
text( 0.0, -3.5, make_label(reg$C,   26), cex = 0.85)
text( 0.0,  3.1, make_label(reg$AB,  20), cex = 0.85)
text(-2.0, -0.9, make_label(reg$AC,  20), cex = 0.85)
text( 2.0, -0.9, make_label(reg$BC,  20), cex = 0.85)
text( 0.0,  0.2, make_label(reg$ABC, 20), cex = 0.9, font = 2)

# circle titles
text(-3.6,  4.7, "Dataset 1\n(Microarray GSE35182)", col = cols["A"], font = 2, cex = 1.3)
text( 3.6,  4.7, "Dataset 2\n(RNA-Seq GSE248521)", col = cols["B"], font = 2, cex = 1.3)
text( 0.0, -5.7, "Dataset 3\n(Microarray GSE19496)", col = cols["C"], font = 2, cex = 1.3)

dev.off()
