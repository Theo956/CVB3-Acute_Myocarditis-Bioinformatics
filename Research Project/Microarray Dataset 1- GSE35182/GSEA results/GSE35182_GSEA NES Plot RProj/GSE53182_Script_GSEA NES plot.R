# Install/load packages -------------------------------------------------------------------------
install.packages("ggplot2")
install.packages("dplyr")
library(ggplot2)
library(dplyr)


# Load iDEP pathway data ----
pathwaydata <- read.csv("GSE53182_Sig_Pathways.csv", stringsAsFactors = FALSE)

names(pathwaydata)

Pathway <- pathwaydata %>%
  rename(
    pathway = matches("Pathways"),   
    NES     = NES,
    pval    = adj.Pval)

top15_Pathway <- head(Pathway[order(Pathway$pval), ], 15)
top15_Pathway

# Bar chart of NES for top15 ----------------------------------------------------------------
ggplot(top15_Pathway, aes(x = NES, y = reorder(pathway, -pval), fill = factor(sign(NES)))) +
  geom_bar(stat = "identity", width = 0.8) +
  labs(title = "GSEA", x = "Normalised Enrichment Score (NES)", y = "Pathway") +
  theme_minimal(base_size = 16) +
  scale_fill_manual(values = c("#B10029", "#0754A2"), guide = "none") +
  scale_y_discrete(labels = function(x) gsub("^KEGG_", "", x)) +
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Save ---------------------------------------------------------------------------------
ggsave("GSE53182_top15_GSEA NES plot.pdf", width = 16, height = 12)
