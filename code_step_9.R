## Load packages

library(clusterProfiler)
library(org.Lmajor.eg.db)
library(forcats)
library(ggplot2)
library(dplyr)

##################### Over Representation Analysis (ORA) #######################

# Lists' creation 

Gene_up_input <- Results_input[Results_input$diffexpressed=="UP",2]
Gene_down_input <- Results_input[Results_input$diffexpressed=="DOWN",2]
universe_input <- Results_input$ID

Gene_up_Mono <- Results_Mono[Results_Mono$diffexpressed=="UP",2]
Gene_down_Mono <- Results_Mono[Results_Mono$diffexpressed=="DOWN",2]
universe_Mono <- Results_Mono$ID

Gene_up_LP <- Results_LP[Results_LP$diffexpressed=="UP",2]
Gene_down_LP <- Results_LP[Results_LP$diffexpressed=="DOWN",2]
universe_LP <- Results_LP$ID

Gene_up_HP <- Results_HP[Results_HP$diffexpressed=="UP",2]
Gene_down_HP <- Results_HP[Results_HP$diffexpressed=="DOWN",2]
universe_HP <- Results_HP$ID



### ORA for inputs - UPREGULATED

ORA_input_up <- enrichGO(gene= Gene_up_input$ID,
                         universe      = universe_input,
                         OrgDb         = org.Lmajor.eg.db,
                         keyType = "GID",
                         ont           = "ALL",
                         pAdjustMethod = "fdr",
                         pvalueCutoff  = 1,
                         qvalueCutoff  = 1,
                         readable      = TRUE)

# Convert enrichment results to a data frame
ORA_input_up_df <- as.data.frame(ORA_input_up)

# Separate and sort top GO terms for each category
Top_BP_up_input <- ORA_input_up_df[ORA_input_up_df$ONTOLOGY == "BP", ][order(ORA_input_up_df[ORA_input_up_df$ONTOLOGY == "BP", ]$pvalue), ][1:15, ]
Top_CC_up_input <- ORA_input_up_df[ORA_input_up_df$ONTOLOGY == "CC", ][order(ORA_input_up_df[ORA_input_up_df$ONTOLOGY == "CC", ]$pvalue), ][1:15, ]
Top_MF_up_input <- ORA_input_up_df[ORA_input_up_df$ONTOLOGY == "MF", ][order(ORA_input_up_df[ORA_input_up_df$ONTOLOGY == "MF", ]$pvalue), ][1:15, ]

# Combine the top terms from each category
Top_ORA_up_input <- rbind(Top_BP_up_input, Top_CC_up_input, Top_MF_up_input)

# Add a column for the GO category
Top_ORA_up_input$Category <- factor(Top_ORA_up_input$ONTOLOGY,
                                    levels = c("BP", "CC", "MF"),
                                    labels= c("Biological process","Cellular component", "Molecular function"))

# Create a custom ordering factor for Description
Top_ORA_up_input$Description <- factor(Top_ORA_up_input$Description, 
                                       levels = Top_ORA_up_input$Description[order(Top_ORA_up_input$Category, -Top_ORA_up_input$pvalue)])

# Plot 

plot <- ggplot(Top_ORA_up_input, aes(x = Description, y = -log(pvalue), fill = Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Molecular function" = "cornflowerblue",
                               "Cellular component" = "darkolivegreen3", 
                               "Biological process" = "lightcoral")) +
  labs(title = "ORA Inputs - upregulated",
       x = "GO term",
       y = "Log(Pvalue)",
       fill = "GO category") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 8),
        legend.position = "bottom")

png("ORA Inputs - upregulated.png", units = "in", width=9, height=10, res=500)
plot
dev.off()

### ORA for inputs - DOWNREGULATED

ORA_input_down <- enrichGO(gene= Gene_down_input$ID,
                           universe      = universe_input,
                           OrgDb         = org.Lmajor.eg.db,
                           keyType = "GID",
                           ont           = "ALL",
                           pAdjustMethod = "fdr",
                           pvalueCutoff  = 1,
                           qvalueCutoff  = 1,
                           readable      = TRUE)

# Convert enrichment results to a data frame
ORA_input_down_df <- as.data.frame(ORA_input_down)

# Separate and sort top GO terms for each category
Top_BP_down_input <- ORA_input_down_df[ORA_input_down_df$ONTOLOGY == "BP", ][order(ORA_input_down_df[ORA_input_down_df$ONTOLOGY == "BP", ]$pvalue), ][1:15, ]
Top_CC_down_input <- ORA_input_down_df[ORA_input_down_df$ONTOLOGY == "CC", ][order(ORA_input_down_df[ORA_input_down_df$ONTOLOGY == "CC", ]$pvalue), ][1:15, ]
Top_MF_down_input <- ORA_input_down_df[ORA_input_down_df$ONTOLOGY == "MF", ][order(ORA_input_down_df[ORA_input_down_df$ONTOLOGY == "MF", ]$pvalue), ][1:15, ]

# Combine the top terms from each category
Top_ORA_down_input <- rbind(Top_BP_down_input, Top_CC_down_input, Top_MF_down_input)

# Add a column for the GO category
Top_ORA_down_input$Category <- factor(Top_ORA_down_input$ONTOLOGY,
                                      levels = c("BP", "CC", "MF"),
                                      labels= c("Biological process","Cellular component", "Molecular function"))

# Create a custom ordering factor for Description
Top_ORA_down_input$Description <- factor(Top_ORA_down_input$Description, 
                                         levels = Top_ORA_down_input$Description[order(Top_ORA_down_input$Category, -Top_ORA_down_input$pvalue)])

# Plot 

plot <- ggplot(Top_ORA_down_input, aes(x = Description, y = -log(pvalue), fill = Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Molecular function" = "cornflowerblue",
                               "Cellular component" = "darkolivegreen3", 
                               "Biological process" = "lightcoral")) +
  labs(title = "ORA Inputs - downregulated",
       x = "GO term",
       y = "Log(Pvalue)",
       fill = "GO category") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 8),
        legend.position = "bottom")

png("ORA Inputs - downregulated.png", units = "in", width=9, height=10, res=500)
plot
dev.off()


### ORA for Monosomes - UPREGULATED

ORA_Mono_up <- enrichGO(gene= Gene_up_Mono$ID,
                        universe      = universe_Mono,
                        OrgDb         = org.Lmajor.eg.db,
                        keyType = "GID",
                        ont           = "ALL",
                        pAdjustMethod = "fdr",
                        pvalueCutoff  = 1,
                        qvalueCutoff  = 1,
                        readable      = TRUE)

# Convert enrichment results to a data frame
ORA_Mono_up_df <- as.data.frame(ORA_Mono_up)

# Separate and sort top GO terms for each category
Top_BP_up_Mono <- ORA_Mono_up_df[ORA_Mono_up_df$ONTOLOGY == "BP", ][order(ORA_Mono_up_df[ORA_Mono_up_df$ONTOLOGY == "BP", ]$pvalue), ][1:15, ]
Top_CC_up_Mono <- ORA_Mono_up_df[ORA_Mono_up_df$ONTOLOGY == "CC", ][order(ORA_Mono_up_df[ORA_Mono_up_df$ONTOLOGY == "CC", ]$pvalue), ][1:15, ]
Top_MF_up_Mono <- ORA_Mono_up_df[ORA_Mono_up_df$ONTOLOGY == "MF", ][order(ORA_Mono_up_df[ORA_Mono_up_df$ONTOLOGY == "MF", ]$pvalue), ][1:15, ]

# Combine the top terms from each category
Top_ORA_up_Mono <- rbind(Top_BP_up_Mono, Top_CC_up_Mono, Top_MF_up_Mono)

# Add a column for the GO category
Top_ORA_up_Mono$Category <- factor(Top_ORA_up_Mono$ONTOLOGY,
                                   levels = c("BP", "CC", "MF"),
                                   labels= c("Biological process","Cellular component", "Molecular function"))

# Create a custom ordering factor for Description
Top_ORA_up_Mono$Description <- factor(Top_ORA_up_Mono$Description, 
                                      levels = Top_ORA_up_Mono$Description[order(Top_ORA_up_Mono$Category, -Top_ORA_up_Mono$pvalue)])

# Plot 

plot <- ggplot(Top_ORA_up_Mono, aes(x = Description, y = -log(pvalue), fill = Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Molecular function" = "cornflowerblue",
                               "Cellular component" = "darkolivegreen3", 
                               "Biological process" = "lightcoral")) +
  labs(title = "ORA Monosomes - upregulated",
       x = "GO term",
       y = "Log(Pvalue)",
       fill = "GO category") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 8),
        legend.position = "bottom")

png("ORA Monosomes - upregulated.png", units = "in", width=10, height=10, res=500)
plot
dev.off()

### ORA for Monosomes - DOWNREGULATED

ORA_Mono_down <- enrichGO(gene= Gene_down_Mono$ID,
                          universe      = universe_Mono,
                          OrgDb         = org.Lmajor.eg.db,
                          keyType = "GID",
                          ont           = "ALL",
                          pAdjustMethod = "fdr",
                          pvalueCutoff  = 1,
                          qvalueCutoff  = 1,
                          readable      = TRUE)

# Convert enrichment results to a data frame
ORA_Mono_down_df <- as.data.frame(ORA_Mono_down)

# Separate and sort top GO terms for each category
Top_BP_down_Mono <- ORA_Mono_down_df[ORA_Mono_down_df$ONTOLOGY == "BP", ][order(ORA_Mono_down_df[ORA_Mono_down_df$ONTOLOGY == "BP", ]$pvalue), ][1:15, ]
Top_CC_down_Mono <- ORA_Mono_down_df[ORA_Mono_down_df$ONTOLOGY == "CC", ][order(ORA_Mono_down_df[ORA_Mono_down_df$ONTOLOGY == "CC", ]$pvalue), ][1:15, ]
Top_MF_down_Mono <- ORA_Mono_down_df[ORA_Mono_down_df$ONTOLOGY == "MF", ][order(ORA_Mono_down_df[ORA_Mono_down_df$ONTOLOGY == "MF", ]$pvalue), ][1:15, ]

# Combine the top terms from each category
Top_ORA_down_Mono <- rbind(Top_BP_down_Mono, Top_CC_down_Mono, Top_MF_down_Mono)

# Add a column for the GO category
Top_ORA_down_Mono$Category <- factor(Top_ORA_down_Mono$ONTOLOGY,
                                     levels = c("BP", "CC", "MF"),
                                     labels= c("Biological process","Cellular component", "Molecular function"))

# Create a custom ordering factor for Description
Top_ORA_down_Mono$Description <- factor(Top_ORA_down_Mono$Description, 
                                        levels = Top_ORA_down_Mono$Description[order(Top_ORA_down_Mono$Category, -Top_ORA_down_Mono$pvalue)])

# Plot 

plot <- ggplot(Top_ORA_down_Mono, aes(x = Description, y = -log(pvalue), fill = Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Molecular function" = "cornflowerblue",
                               "Cellular component" = "darkolivegreen3", 
                               "Biological process" = "lightcoral")) +
  labs(title = "ORA Monosomes - downregulated",
       x = "GO term",
       y = "Log(Pvalue)",
       fill = "GO category") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 8),
        legend.position = "bottom")

png("ORA Monosomes - downregulated.png", units = "in", width=10, height=10, res=500)
plot
dev.off()


### ORA for Light polysomes - UPREGULATED

ORA_LP_up <- enrichGO(gene= Gene_up_LP$ID,
                      universe      = universe_LP,
                      OrgDb         = org.Lmajor.eg.db,
                      keyType = "GID",
                      ont           = "ALL",
                      pAdjustMethod = "fdr",
                      pvalueCutoff  = 1,
                      qvalueCutoff  = 1,
                      readable      = TRUE)

# Convert enrichment results to a data frame
ORA_LP_up_df <- as.data.frame(ORA_LP_up)

# Separate and sort top GO terms for each category
Top_BP_up_LP <- ORA_LP_up_df[ORA_LP_up_df$ONTOLOGY == "BP", ][order(ORA_LP_up_df[ORA_LP_up_df$ONTOLOGY == "BP", ]$pvalue), ][1:15, ]
Top_CC_up_LP <- ORA_LP_up_df[ORA_LP_up_df$ONTOLOGY == "CC", ][order(ORA_LP_up_df[ORA_LP_up_df$ONTOLOGY == "CC", ]$pvalue), ][1:15, ]
Top_MF_up_LP <- ORA_LP_up_df[ORA_LP_up_df$ONTOLOGY == "MF", ][order(ORA_LP_up_df[ORA_LP_up_df$ONTOLOGY == "MF", ]$pvalue), ][1:15, ]

# Combine the top terms from each category
Top_ORA_up_LP <- rbind(Top_BP_up_LP, Top_CC_up_LP, Top_MF_up_LP)

# Add a column for the GO category
Top_ORA_up_LP$Category <- factor(Top_ORA_up_LP$ONTOLOGY,
                                 levels = c("BP", "CC", "MF"),
                                 labels= c("Biological process","Cellular component", "Molecular function"))

# Create a custom ordering factor for Description
Top_ORA_up_LP$Description <- factor(Top_ORA_up_LP$Description, 
                                    levels = Top_ORA_up_LP$Description[order(Top_ORA_up_LP$Category, -Top_ORA_up_LP$pvalue)])

# Plot 

plot <- ggplot(Top_ORA_up_LP, aes(x = Description, y = -log(pvalue), fill = Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Molecular function" = "cornflowerblue",
                               "Cellular component" = "darkolivegreen3", 
                               "Biological process" = "lightcoral")) +
  labs(title = "ORA Light polysomes - upregulated",
       x = "GO term",
       y = "Log(Pvalue)",
       fill = "GO category") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 8),
        legend.position = "bottom")

png("ORA Light polysomes - upregulated.png", units = "in", width=9, height=10, res=500)
plot
dev.off()

### ORA for Light polysomes - DOWNREGULATED

ORA_LP_down <- enrichGO(gene= Gene_down_LP$ID,
                        universe      = universe_LP,
                        OrgDb         = org.Lmajor.eg.db,
                        keyType = "GID",
                        ont           = "ALL",
                        pAdjustMethod = "fdr",
                        pvalueCutoff  = 1,
                        qvalueCutoff  = 1,
                        readable      = TRUE)

# Convert enrichment results to a data frame
ORA_LP_down_df <- as.data.frame(ORA_LP_down)

# Separate and sort top GO terms for each category
Top_BP_down_LP <- ORA_LP_down_df[ORA_LP_down_df$ONTOLOGY == "BP", ][order(ORA_LP_down_df[ORA_LP_down_df$ONTOLOGY == "BP", ]$pvalue), ][1:15, ]
Top_CC_down_LP <- ORA_LP_down_df[ORA_LP_down_df$ONTOLOGY == "CC", ][order(ORA_LP_down_df[ORA_LP_down_df$ONTOLOGY == "CC", ]$pvalue), ][1:15, ]
Top_MF_down_LP <- ORA_LP_down_df[ORA_LP_down_df$ONTOLOGY == "MF", ][order(ORA_LP_down_df[ORA_LP_down_df$ONTOLOGY == "MF", ]$pvalue), ][1:15, ]

# Combine the top terms from each category
Top_ORA_down_LP <- rbind(Top_BP_down_LP, Top_CC_down_LP, Top_MF_down_LP)

# Add a column for the GO category
Top_ORA_down_LP$Category <- factor(Top_ORA_down_LP$ONTOLOGY,
                                   levels = c("BP", "CC", "MF"),
                                   labels= c("Biological process","Cellular component", "Molecular function"))

# Create a custom ordering factor for Description
Top_ORA_down_LP$Description <- factor(Top_ORA_down_LP$Description, 
                                      levels = Top_ORA_down_LP$Description[order(Top_ORA_down_LP$Category, -Top_ORA_down_LP$pvalue)])

# Plot 

plot <- ggplot(Top_ORA_down_LP, aes(x = Description, y = -log(pvalue), fill = Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Molecular function" = "cornflowerblue",
                               "Cellular component" = "darkolivegreen3", 
                               "Biological process" = "lightcoral")) +
  labs(title = "ORA Light polysomes - downregulated",
       x = "GO term",
       y = "Log(Pvalue)",
       fill = "GO category") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 8),
        legend.position = "bottom")

png("ORA Light polysomes - downregulated.png", units = "in", width=9, height=10, res=500)
plot
dev.off()


### ORA for Heavy polysomes - UPREGULATED

ORA_HP_up <- enrichGO(gene= Gene_up_HP$ID,
                      universe      = universe_HP,
                      OrgDb         = org.Lmajor.eg.db,
                      keyType = "GID",
                      ont           = "ALL",
                      pAdjustMethod = "fdr",
                      pvalueCutoff  = 1,
                      qvalueCutoff  = 1,
                      readable      = TRUE)

# Convert enrichment results to a data frame
ORA_HP_up_df <- as.data.frame(ORA_HP_up)

# Separate and sort top GO terms for each category
Top_BP_up_HP <- ORA_HP_up_df[ORA_HP_up_df$ONTOLOGY == "BP", ][order(ORA_HP_up_df[ORA_HP_up_df$ONTOLOGY == "BP", ]$pvalue), ][1:15, ]
Top_CC_up_HP <- ORA_HP_up_df[ORA_HP_up_df$ONTOLOGY == "CC", ][order(ORA_HP_up_df[ORA_HP_up_df$ONTOLOGY == "CC", ]$pvalue), ][1:15, ]
Top_MF_up_HP <- ORA_HP_up_df[ORA_HP_up_df$ONTOLOGY == "MF", ][order(ORA_HP_up_df[ORA_HP_up_df$ONTOLOGY == "MF", ]$pvalue), ][1:15, ]

# Combine the top terms from each category
Top_ORA_up_HP <- rbind(Top_BP_up_HP, Top_CC_up_HP, Top_MF_up_HP)

# Add a column for the GO category
Top_ORA_up_HP$Category <- factor(Top_ORA_up_HP$ONTOLOGY,
                                 levels = c("BP", "CC", "MF"),
                                 labels= c("Biological process","Cellular component", "Molecular function"))

# Create a custom ordering factor for Description
Top_ORA_up_HP$Description <- factor(Top_ORA_up_HP$Description, 
                                    levels = Top_ORA_up_HP$Description[order(Top_ORA_up_HP$Category, -Top_ORA_up_HP$pvalue)])

# Plot 

plot <- ggplot(Top_ORA_up_HP, aes(x = Description, y = -log(pvalue), fill = Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Molecular function" = "cornflowerblue",
                               "Cellular component" = "darkolivegreen3", 
                               "Biological process" = "lightcoral")) +
  labs(title = "ORA Heavy polysomes - upregulated",
       x = "GO term",
       y = "Log(Pvalue)",
       fill = "GO category") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 8),
        legend.position = "bottom")

png("ORA Heavy polysomes - upregulated.png", units = "in", width=9, height=10, res=500)
plot
dev.off()

### ORA for Heavy polysomes - DOWNREGULATED

ORA_HP_down <- enrichGO(gene= Gene_down_HP$ID,
                        universe      = universe_HP,
                        OrgDb         = org.Lmajor.eg.db,
                        keyType = "GID",
                        ont           = "ALL",
                        pAdjustMethod = "fdr",
                        pvalueCutoff  = 1,
                        qvalueCutoff  = 1,
                        readable      = TRUE)

# Convert enrichment results to a data frame
ORA_HP_down_df <- as.data.frame(ORA_HP_down)

# Separate and sort top GO terms for each category
Top_BP_down_HP <- ORA_HP_down_df[ORA_HP_down_df$ONTOLOGY == "BP", ][order(ORA_HP_down_df[ORA_HP_down_df$ONTOLOGY == "BP", ]$pvalue), ][1:15, ]
Top_CC_down_HP <- ORA_HP_down_df[ORA_HP_down_df$ONTOLOGY == "CC", ][order(ORA_HP_down_df[ORA_HP_down_df$ONTOLOGY == "CC", ]$pvalue), ][1:15, ]
Top_MF_down_HP <- ORA_HP_down_df[ORA_HP_down_df$ONTOLOGY == "MF", ][order(ORA_HP_down_df[ORA_HP_down_df$ONTOLOGY == "MF", ]$pvalue), ][1:15, ]

# Combine the top terms from each category
Top_ORA_down_HP <- rbind(Top_BP_down_HP, Top_CC_down_HP, Top_MF_down_HP)

# Add a column for the GO category
Top_ORA_down_HP$Category <- factor(Top_ORA_down_HP$ONTOLOGY,
                                   levels = c("BP", "CC", "MF"),
                                   labels= c("Biological process","Cellular component", "Molecular function"))

# Create a custom ordering factor for Description
Top_ORA_down_HP$Description <- factor(Top_ORA_down_HP$Description, 
                                      levels = Top_ORA_down_HP$Description[order(Top_ORA_down_HP$Category, -Top_ORA_down_HP$pvalue)])

# Plot 

plot <- ggplot(Top_ORA_down_HP, aes(x = Description, y = -log(pvalue), fill = Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Molecular function" = "cornflowerblue",
                               "Cellular component" = "darkolivegreen3", 
                               "Biological process" = "lightcoral")) +
  labs(title = "ORA Heavy polysomes - downregulated",
       x = "GO term",
       y = "Log(Pvalue)",
       fill = "GO category") +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 8),
        legend.position = "bottom")

png("ORA Heavy polysomes - downregulated.png", units = "in", width=9, height=10, res=500)
plot
dev.off()


################## Gene set enrichment analysis (GSEA) ########################

# Lists' creation

GSEA_input_gene_list <- Results_input$logFC
names(GSEA_input_gene_list) <- Results_input$ID
GSEA_input_gene_list <- na.omit(GSEA_input_gene_list)
GSEA_input_gene_list <- sort(GSEA_input_gene_list, decreasing = T)

GSEA_Mono_gene_list <- Results_Mono$logFC
names(GSEA_Mono_gene_list) <- Results_Mono$ID
GSEA_Mono_gene_list <- na.omit(GSEA_Mono_gene_list)
GSEA_Mono_gene_list <- sort(GSEA_Mono_gene_list, decreasing = T)

GSEA_LP_gene_list <- Results_LP$logFC
names(GSEA_LP_gene_list) <- Results_LP$ID
GSEA_LP_gene_list <- na.omit(GSEA_LP_gene_list)
GSEA_LP_gene_list <- sort(GSEA_LP_gene_list, decreasing = T)

GSEA_HP_gene_list <- Results_HP$logFC
names(GSEA_HP_gene_list) <- Results_HP$ID
GSEA_HP_gene_list <- na.omit(GSEA_HP_gene_list)
GSEA_HP_gene_list <- sort(GSEA_HP_gene_list, decreasing = T)


# Function to get top 10 upregulated and downregulated terms for a specific category
get_top_terms <- function(data, category) {
  upregulated_terms <- data %>%
    filter(ONTOLOGY == category, NES > 0) %>%
    slice_max(order_by = NES, n = 10)
  
  downregulated_terms <- data %>%
    filter(ONTOLOGY == category, NES < 0) %>%
    slice_min(order_by = NES, n = 10)
  
  combined_terms <- bind_rows(upregulated_terms, downregulated_terms) %>%
    mutate(term = fct_reorder(Description, NES),
           category = category)
  
  return(combined_terms)
}

### Inputs

GSEA_input <- gseGO(geneList = GSEA_input_gene_list,
                    ont = "ALL",
                    keyType = "GID",
                    exponent=1,
                    minGSSize = 10,
                    maxGSSize = 200,
                    eps = 1e-10,
                    pvalueCutoff = 1,
                    OrgDb = "org.Lmajor.eg.db",
                    pAdjustMethod = "fdr")

#Convert GSEA output into a readable dataframe
GSEA_input_df <- as.data.frame(GSEA_input)

# Get the top terms for inputs
BP_terms_input <- get_top_terms(GSEA_input@result, "BP")
CC_terms_input <- get_top_terms(GSEA_input@result, "CC")
MF_terms_input <- get_top_terms(GSEA_input@result, "MF")

# Combine all terms
All_terms_input <- bind_rows(BP_terms_input, CC_terms_input, MF_terms_input)

# Plot
GSEA_input_barplot <- ggplot(All_terms_input, aes(term, NES, fill = NES)) +
  geom_col() +
  coord_flip() +
  theme_test() +
  scale_fill_gradient(high = "firebrick2", low = "royalblue1") +
  geom_text(aes(label = formatC(pvalue, format = "e", digits = 2)),
            hjust = ifelse(All_terms_input$NES > 0, -0.1, 1.1), color = "white") +
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_text(size = 8)) +
  geom_hline(yintercept = 0) +
  ggtitle("GSEA - Inputs")+
  ylab("NES") +
  xlab("GO term") +
  facet_wrap(~category, scales = "free_y", ncol = 1)
png("GSEA - Inputs.png", units = "in", width=9, height=10, res=500)
GSEA_input_barplot
dev.off()

### Monosomes

GSEA_Mono <- gseGO(geneList = GSEA_Mono_gene_list,
                   ont = "ALL",
                   keyType = "GID",
                   exponent=1,
                   minGSSize = 10,
                   maxGSSize = 200,
                   eps = 1e-10,
                   pvalueCutoff = 1,
                   OrgDb = "org.Lmajor.eg.db",
                   pAdjustMethod = "fdr")

#Convert GSEA output into a readable dataframe
GSEA_Mono_df <- as.data.frame(GSEA_Mono)

# Get the top terms for Monosomes
BP_terms_Mono <- get_top_terms(GSEA_Mono@result, "BP")
CC_terms_Mono <- get_top_terms(GSEA_Mono@result, "CC")
MF_terms_Mono <- get_top_terms(GSEA_Mono@result, "MF")

# Combine all terms
All_terms_Mono <- bind_rows(BP_terms_Mono, CC_terms_Mono, MF_terms_Mono)

# Plot
GSEA_Mono_barplot <- ggplot(All_terms_Mono, aes(term, NES, fill = NES)) +
  geom_col() +
  coord_flip() +
  theme_test() +
  scale_fill_gradient(high = "firebrick2", low = "royalblue1") +
  geom_text(aes(label = formatC(pvalue, format = "e", digits = 2)),
            hjust = ifelse(All_terms_Mono$NES > 0, -0.1, 1.1), color = "white") +
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_text(size = 8)) +
  geom_hline(yintercept = 0) +
  ggtitle("GSEA - Monosomes")+
  ylab("NES") +
  xlab("GO term") +
  facet_wrap(~category, scales = "free_y", ncol = 1)
png("GSEA - Monosomes.png", units = "in", width=9, height=10, res=500)
GSEA_Mono_barplot
dev.off()

### Light polysomes

GSEA_LP <- gseGO(geneList = GSEA_LP_gene_list,
                 ont = "ALL",
                 keyType = "GID",
                 exponent=1,
                 minGSSize = 10,
                 maxGSSize = 200,
                 eps = 1e-10,
                 pvalueCutoff = 1,
                 OrgDb = "org.Lmajor.eg.db",
                 pAdjustMethod = "fdr")

#Convert GSEA output into a readable dataframe
GSEA_LP_df <- as.data.frame(GSEA_LP)

# Get the top terms for Light polysomes
BP_terms_LP <- get_top_terms(GSEA_LP@result, "BP")
CC_terms_LP <- get_top_terms(GSEA_LP@result, "CC")
MF_terms_LP <- get_top_terms(GSEA_LP@result, "MF")

# Combine all terms
All_terms_LP <- bind_rows(BP_terms_LP, CC_terms_LP, MF_terms_LP)

# Plot
GSEA_LP_barplot <- ggplot(All_terms_LP, aes(term, NES, fill = NES)) +
  geom_col() +
  coord_flip() +
  theme_test() +
  scale_fill_gradient(high = "firebrick2", low = "royalblue1") +
  geom_text(aes(label = formatC(pvalue, format = "e", digits = 2)),
            hjust = ifelse(All_terms_LP$NES > 0, -0.1, 1.1), color = "white") +
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_text(size = 8)) +
  geom_hline(yintercept = 0) +
  ggtitle("GSEA - Light polysomes")+
  ylab("NES") +
  xlab("GO term") +
  facet_wrap(~category, scales = "free_y", ncol = 1)
png("GSEA - Light polysomes.png", units = "in", width=9, height=10, res=500)
GSEA_LP_barplot
dev.off()


### Heavy polysomes

GSEA_HP <- gseGO(geneList = GSEA_HP_gene_list,
                 ont = "ALL",
                 keyType = "GID",
                 exponent=1,
                 minGSSize = 10,
                 maxGSSize = 200,
                 eps = 1e-10,
                 pvalueCutoff = 1,
                 OrgDb = "org.Lmajor.eg.db",
                 pAdjustMethod = "fdr")

#Convert GSEA output into a readable dataframe
GSEA_HP_df <- as.data.frame(GSEA_HP)

# Get the top terms for Heavy polysomes
BP_terms_HP <- get_top_terms(GSEA_HP@result, "BP")
CC_terms_HP <- get_top_terms(GSEA_HP@result, "CC")
MF_terms_HP <- get_top_terms(GSEA_HP@result, "MF")

# Combine all terms
All_terms_HP <- bind_rows(BP_terms_HP, CC_terms_HP, MF_terms_HP)

# Plot
GSEA_HP_barplot <- ggplot(All_terms_HP, aes(term, NES, fill = NES)) +
  geom_col() +
  coord_flip() +
  theme_test() +
  scale_fill_gradient(high = "firebrick2", low = "royalblue1") +
  geom_text(aes(label = formatC(pvalue, format = "e", digits = 2)),
            hjust = ifelse(All_terms_HP$NES > 0, -0.1, 1.1), color = "white") +
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_text(size = 8)) +
  geom_hline(yintercept = 0) +
  ggtitle("GSEA - Heavy polysomes")+
  ylab("NES") +
  xlab("GO term") +
  facet_wrap(~category, scales = "free_y", ncol = 1)
png("GSEA - Heavy polysomes.png", units = "in", width=9, height=10, res=500)
GSEA_HP_barplot
dev.off()
