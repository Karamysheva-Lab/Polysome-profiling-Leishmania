library(edgeR)
library(ggfortify)
library(ggplot2)
library(corrplot)
library(ggrepel)


## Import raw counts and sample table

raw_counts <- read.delim("/Volumes/RAID5_Ponomarev/Users/chrisbu/raw_counts.txt",
                         comment.char="#")
sample_table <- read.delim("/Volumes/RAID5_Ponomarev/Users/chrisbu/sample_table.txt", header = F)

## Pre-processing raw counts and sample table

raw_counts$Geneid <- sub("^exon_", "", raw_counts$Geneid)
raw_counts$Geneid <- sub(":.*", "", raw_counts$Geneid)
raw_counts <- raw_counts[,c(-2,-3,-4,-5,-6)]

table(duplicated(raw_counts$Geneid))
colnames(raw_counts)[2:25] <-  1:24
sample_table$sample <- 1:24
colnames(sample_table)[1] <- "sample_name"
sample_table$fraction <- c(rep("input",6),rep("Mono",6),rep("LP",6),rep("HP",6))
sample_table$group <- c(rep(c("control","control","control","treat","treat","treat"),
                              4))
raw_counts$Geneid <- as.factor(raw_counts$Geneid)
raw_counts <- aggregate(. ~Geneid, data = raw_counts, FUN = sum)  
rownames(raw_counts) <- raw_counts$Geneid
raw_counts <- raw_counts[,-1]

## Create a edgeR object with all samples

DGE <- DGEList(counts = raw_counts, remove.zeros = TRUE,
               group = sample_table$group)
DGE$samples

#################  PCA  and correlation with all samples ####################

# get normalized counts
CPM <- cpm(DGE,normalized.lib.sizes = TRUE, log = TRUE)
#PCA
PCA <- prcomp(t(CPM))
#plot
png("PCA_plot_1.png", units = "in", width=7, height=6, res=500)
autoplot(PCA,x=1, y=2, colour = "group", shape = "fraction",size=3,data = sample_table) + 
  geom_hline(yintercept=0, col="gray") + geom_vline(xintercept=0,col="gray") + theme_minimal() + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) + 
  guides(colour=guide_legend(title="Group"),shape=guide_legend(title="Time"))
dev.off()

# correlation
Cor_matrix <- cor(CPM)
#plot
png("Cor_plot_1.png", units = "in", width=7, height=6, res=500)
corrplot(Cor_matrix, method = 'color',is.corr = FALSE, col.lim = c(0.7, 1), col = COL1("Blues"))
dev.off()

#################  exploratory analysis with all samples ####################

# barplot for all samples
png("barplot_raw.png", units = "in", width=8, height=6, res=500)
barplot(colSums(raw_counts), las=1,cex.names=0.8)
dev.off()

# boxplot by fraction and group for all samples
boxplot(colSums(raw_counts) ~ sample_table$fraction)
boxplot(colSums(raw_counts) ~ sample_table$group)

MinVals <- apply(raw_counts, 1, min)
sum(MinVals == 0)
Exp <- raw_counts[MinVals > 0, ]
Exp_log2 <- as.matrix(log2(Exp))

# boxplot of raw counts by sample
png("boxplot_raw.png", units = "in", width=8, height=6, res=500)
boxplot(Exp_log2, ylab="log2 counts", main="Raw Counts", las=1, cex.axis=0.8)
dev.off()

DGE <- calcNormFactors(DGE)
DGE$samples
CPM_size <- cpm(DGE,normalized.lib.sizes = TRUE, log = TRUE)

MinVals <- apply(CPM_size, 1, min)
sum(MinVals == 0)
Exp <- CPM_size[MinVals > 0, ]
Exp_log2 <- as.matrix(Exp)

# boxplot of raw counts by sample
png("boxplot_norm.png", units = "in", width=8, height=6, res=500)
boxplot(Exp_log2, ylab="log2 counts", main="Normalized Counts", las=1, cex.axis=0.8)
dev.off()

################# split counts and sample table by fraction  #################

raw_counts_input <- raw_counts[,sample_table$fraction=="input"]
raw_counts_Mono <- raw_counts[,sample_table$fraction=="Mono"]
raw_counts_LP <- raw_counts[,sample_table$fraction=="LP"]
raw_counts_HP <- raw_counts[,sample_table$fraction=="HP"]

sample_table_input <- sample_table[sample_table$fraction=="input",]
sample_table_Mono <- sample_table[sample_table$fraction=="Mono",]
sample_table_LP <- sample_table[sample_table$fraction=="LP",]
sample_table_HP <- sample_table[sample_table$fraction=="HP",]


###############################################################
################## Analysis for INPUT  ########################
###############################################################


# edgeR analysis 
DGE_input <- DGEList(counts = raw_counts_input, remove.zeros = TRUE,
                     group = sample_table_input$group)
design_input <- model.matrix(~0 + sample_table_input$group)
rownames(design_input) <- colnames(DGE_input)
colnames(design_input) <- c("control","treat")
keep <- filterByExpr(DGE_input, design_input, min.count = 5)
DGE_input <- DGE_input[keep, ] 
DGE_input <- calcNormFactors(DGE_input)
DGE_input$samples
DGE_input <- estimateGLMCommonDisp(DGE_input, design_input, verbose = TRUE) 
DGE_input <- estimateGLMTrendedDisp(DGE_input, design_input)
DGE_input <- estimateGLMTagwiseDisp(DGE_input, design_input)
fit_input <- glmQLFit(DGE_input, design_input)
input_cont <- makeContrasts(treat-control, levels = design_input)
QLtest_input <- glmQLFTest(fit_input, contrast=input_cont)
Results_input <- as.data.frame(topTags(QLtest_input, n = dim(DGE_input)[1]))

# Volcano plot 
Results_input$diffexpressed <- "NO"
Results_input$diffexpressed[Results_input$logFC > 0.5849625 & Results_input$FDR < 0.15] <- "UP"
Results_input$diffexpressed[Results_input$logFC < -0.5849625 & Results_input$FDR < 0.15] <- "DOWN"
Results_input <- Results_input %>% 
  group_by(diffexpressed) %>%
  arrange(ifelse(diffexpressed == "DOWN", logFC, -logFC))%>%
  mutate(Top = ifelse(row_number() <= 10, Symbol,""))%>% # Change the top number and the name of the column with the names or symbols
  mutate(Top = ifelse(diffexpressed == "NO", "", Top)) # With this code it will erase the top names or symbols created in NO category
row.names(Results_input) <- Results_input$ID

plot <- ggplot(data = Results_input, aes(x = logFC, y = -log10(PValue), col = diffexpressed, label =Top)) +
  geom_point() +
  theme_test() +
  geom_text_repel(size=4, max.overlaps = 20, nudge_x = .1,
                  box.padding = 0.5,
                  point.padding = 0.01,
                  nudge_y = .3,
                  min.segment.length = 0,
                  segment.size= 0.5,
                  show.legend = FALSE) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "gray", "UP" = "red"),
                     breaks = c("DOWN", "UP"),
                     labels = c("DOWN", "UP")) +
  geom_vline(xintercept = -1, col = "gray30", linetype = "dashed") +
  geom_vline(xintercept = 1, col = "gray30", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.003847380), col = "gray30", linetype = "dashed") +
  theme(legend.title = element_blank(),axis.text=element_text(size=13), 
        axis.title=element_text(size=15),text =element_text(size=13))
png("Volcano_input.png", units = "in", width=8, height=6, res=500)
plot
dev.off()


CPM_input <- cpm(DGE_input,normalized.lib.sizes = TRUE, log = TRUE)

MinVals <- apply(CPM_input, 1, min)
sum(MinVals == 0)
Exp <- CPM_input[MinVals > 0, ]
Exp_log2 <- as.matrix(Exp)

png("boxplot_input.png", units = "in", width=8, height=6, res=500)
boxplot(Exp_log2, ylab="log2 counts", main="Normalized Counts Input", las=1, cex.axis=0.8)
dev.off()

PCA <- prcomp(t(CPM_input))
png("PCA_input.png", units = "in", width=7, height=6, res=500)
autoplot(PCA,x=1, y=2, colour = "group", shape = "fraction",size=3,data = sample_table_input) + 
  geom_hline(yintercept=0, col="gray") + geom_vline(xintercept=0,col="gray") + theme_minimal() + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) + 
  guides(colour=guide_legend(title="Group"),shape=guide_legend(title="Time"))
dev.off()

Cor_matrix <- cor(CPM_input)
png("Cor_plot_input.png", units = "in", width=7, height=6, res=500)
corrplot(Cor_matrix, method = 'color',is.corr = FALSE, col.lim = c(0.9, 1), col = COL1("Blues"))
dev.off()

###############################################################
################## Analysis for MONO  #########################
###############################################################

# edgeR analysis 
DGE_Mono <- DGEList(counts = raw_counts_Mono, remove.zeros = TRUE,
                     group = sample_table_Mono$group)
design_Mono <- model.matrix(~0 + sample_table_Mono$group)
rownames(design_Mono) <- colnames(DGE_Mono)
colnames(design_Mono) <- c("control","treat")
keep <- filterByExpr(DGE_Mono, design_Mono, min.count = 5)
DGE_Mono <- DGE_Mono[keep, ] 
DGE_Mono <- calcNormFactors(DGE_Mono)
DGE_Mono$samples
DGE_Mono <- estimateGLMCommonDisp(DGE_Mono, design_Mono, verbose = TRUE) 
DGE_Mono <- estimateGLMTrendedDisp(DGE_Mono, design_Mono)
DGE_Mono <- estimateGLMTagwiseDisp(DGE_Mono, design_Mono)
fit_Mono <- glmQLFit(DGE_Mono, design_Mono)
Mono_cont <- makeContrasts(treat-control, levels = design_Mono)
QLtest_Mono <- glmQLFTest(fit_Mono, contrast=Mono_cont)
Results_Mono <- as.data.frame(topTags(QLtest_Mono, n = dim(DGE_Mono)[1]))

Results_Mono$diffexpressed <- "NO"
Results_Mono$diffexpressed[Results_Mono$logFC > 0.5849625 & Results_Mono$FDR < 0.15] <- "UP"
Results_Mono$diffexpressed[Results_Mono$logFC < -0.5849625 & Results_Mono$FDR < 0.15] <- "DOWN"
Results_Mono <- Results_Mono %>% 
  group_by(diffexpressed) %>%
  arrange(ifelse(diffexpressed == "DOWN", logFC, -logFC))%>%
  mutate(Top = ifelse(row_number() <= 10, Symbol,""))%>% 
  mutate(Top = ifelse(diffexpressed == "NO", "", Top)) 
row.names(Results_Mono) <- Results_Mono$ID

plot <- ggplot(data = Results_Mono, aes(x = logFC, y = -log10(PValue), col = diffexpressed, label =Top)) +
  geom_point() +
  theme_test() +
  geom_text_repel(size=4, max.overlaps = 20, nudge_x = .1,
                  box.padding = 0.5,
                  point.padding = 0.01,
                  nudge_y = .3,
                  min.segment.length = 0,
                  segment.size= 0.5,
                  show.legend = FALSE) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "gray", "UP" = "red"),
                     breaks = c("DOWN", "UP"),
                     labels = c("DOWN", "UP")) +
  geom_vline(xintercept = -1, col = "gray30", linetype = "dashed") +
  geom_vline(xintercept = 1, col = "gray30", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.0006068185), col = "gray30", linetype = "dashed") +
  theme(legend.title = element_blank(),axis.text=element_text(size=13), 
        axis.title=element_text(size=15),text =element_text(size=13))
png("Volcano_Mono.png", units = "in", width=8, height=6, res=500)
plot
dev.off()

CPM_Mono <- cpm(DGE_Mono,normalized.lib.sizes = TRUE, log = TRUE)

MinVals <- apply(CPM_Mono, 1, min)
sum(MinVals == 0)
Exp <- CPM_Mono[MinVals > 0, ]
Exp_log2 <- as.matrix(Exp)

png("boxplot_Mono.png", units = "in", width=8, height=6, res=500)
boxplot(Exp_log2, ylab="log2 counts", main="Normalized Counts Mono", las=1, cex.axis=0.8)
dev.off()

PCA <- prcomp(t(CPM_Mono))
png("PCA_Mono.png", units = "in", width=7, height=6, res=500)
autoplot(PCA,label = T,x=1, y=2, colour = "group", shape = "fraction",size=1,data = sample_table_Mono) + 
  geom_hline(yintercept=0, col="gray") + geom_vline(xintercept=0,col="gray") + theme_minimal() + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) + 
  guides(colour=guide_legend(title="Group"),shape=guide_legend(title="Time"))
dev.off()

Cor_matrix <- cor(CPM_Mono)
png("Cor_plot_Mono.png", units = "in", width=7, height=6, res=500)
corrplot(Cor_matrix, method = 'color',is.corr = FALSE, col.lim = c(0.7, 1), col = COL1("Blues"))
dev.off()

###############################################################
################## Analysis for LP  ###########################
###############################################################


# edgeR analysis 
DGE_LP <- DGEList(counts = raw_counts_LP, remove.zeros = TRUE,
                     group = sample_table_LP$group)
design_LP <- model.matrix(~0 + sample_table_LP$group)
rownames(design_LP) <- colnames(DGE_LP)
colnames(design_LP) <- c("control","treat")
keep <- filterByExpr(DGE_LP, design_LP, min.count = 5)
DGE_LP <- DGE_LP[keep, ] 
DGE_LP <- calcNormFactors(DGE_LP)
DGE_LP$samples
DGE_LP <- estimateGLMCommonDisp(DGE_LP, design_LP, verbose = TRUE) 
DGE_LP <- estimateGLMTrendedDisp(DGE_LP, design_LP)
DGE_LP <- estimateGLMTagwiseDisp(DGE_LP, design_LP)
fit_LP <- glmQLFit(DGE_LP, design_LP)
LP_cont <- makeContrasts(treat-control, levels = design_LP)
QLtest_LP <- glmQLFTest(fit_LP, contrast=LP_cont)
Results_LP <- as.data.frame(topTags(QLtest_LP, n = dim(DGE_LP)[1]))

# Volcano plot

Results_LP$diffexpressed <- "NO"
Results_LP$diffexpressed[Results_LP$logFC > 0.5849625 & Results_LP$FDR < 0.15] <- "UP"
Results_LP$diffexpressed[Results_LP$logFC < -0.5849625 & Results_LP$FDR < 0.15] <- "DOWN"
Results_LP <- Results_LP %>% 
  group_by(diffexpressed) %>%
  arrange(ifelse(diffexpressed == "DOWN", logFC, -logFC))%>%
  mutate(Top = ifelse(row_number() <= 10, Symbol,""))%>% 
  mutate(Top = ifelse(diffexpressed == "NO", "", Top)) 
row.names(Results_LP) <- Results_LP$ID

plot <- ggplot(data = Results_LP, aes(x = logFC, y = -log10(PValue), col = diffexpressed, label =Top)) +
  geom_point() +
  theme_test() +
  geom_text_repel(size=4, max.overlaps = 20, nudge_x = .1,
                  box.padding = 0.5,
                  point.padding = 0.01,
                  nudge_y = .3,
                  min.segment.length = 0,
                  segment.size= 0.5,
                  show.legend = FALSE) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "gray", "UP" = "red"),
                     breaks = c("DOWN", "UP"),
                     labels = c("DOWN", "UP")) +
  geom_vline(xintercept = -1, col = "gray30", linetype = "dashed") +
  geom_vline(xintercept = 1, col = "gray30", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.003829267), col = "gray30", linetype = "dashed") +
  theme(legend.title = element_blank(),axis.text=element_text(size=13), 
        axis.title=element_text(size=15),text =element_text(size=13))
png("Volcano_LP.png", units = "in", width=8, height=6, res=500)
plot
dev.off()

CPM_LP <- cpm(DGE_LP,normalized.lib.sizes = TRUE, log = TRUE)

MinVals <- apply(CPM_LP, 1, min)
sum(MinVals == 0)
Exp <- CPM_LP[MinVals > 0, ]
Exp_log2 <- as.matrix(Exp)

png("boxplot_LP.png", units = "in", width=8, height=6, res=500)
boxplot(Exp_log2, ylab="log2 counts", main="Normalized Counts LP", las=1, cex.axis=0.8)
dev.off()

PCA <- prcomp(t(CPM_LP))
png("PCA_LP.png", units = "in", width=7, height=6, res=500)
autoplot(PCA,x=1, y=2, colour = "group", shape = "fraction",size=3,data = sample_table_LP) + 
  geom_hline(yintercept=0, col="gray") + geom_vline(xintercept=0,col="gray") + theme_minimal() + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) + 
  guides(colour=guide_legend(title="Group"),shape=guide_legend(title="Time"))
dev.off()

Cor_matrix <- cor(CPM_LP)
png("Cor_plot_LP.png", units = "in", width=7, height=6, res=500)
corrplot(Cor_matrix, method = 'color',is.corr = FALSE, col.lim = c(0.8, 1), col = COL1("Blues"))
dev.off()


###############################################################
################## Analysis for HP  ###########################
###############################################################

# edgeR analysis 
DGE_HP <- DGEList(counts = raw_counts_HP, remove.zeros = TRUE,
                     group = sample_table_HP$group)
design_HP <- model.matrix(~0 + sample_table_HP$group)
rownames(design_HP) <- colnames(DGE_HP)
colnames(design_HP) <- c("control","treat")
keep <- filterByExpr(DGE_HP, design_HP, min.count = 5)
DGE_HP <- DGE_HP[keep, ] 
DGE_HP <- calcNormFactors(DGE_HP)
DGE_HP$samples
DGE_HP <- estimateGLMCommonDisp(DGE_HP, design_HP, verbose = TRUE) 
DGE_HP <- estimateGLMTrendedDisp(DGE_HP, design_HP)
DGE_HP <- estimateGLMTagwiseDisp(DGE_HP, design_HP)
fit_HP <- glmQLFit(DGE_HP, design_HP)
HP_cont <- makeContrasts(treat-control, levels = design_HP)
QLtest_HP <- glmQLFTest(fit_HP, contrast=HP_cont)
Results_HP <- as.data.frame(topTags(QLtest_HP, n = dim(DGE_HP)[1]))

# Volcano plot
Results_HP$diffexpressed <- "NO"
Results_HP$diffexpressed[Results_HP$logFC > 0.5849625 & Results_HP$FDR < 0.15] <- "UP"
Results_HP$diffexpressed[Results_HP$logFC < -0.5849625 & Results_HP$FDR < 0.15] <- "DOWN"
Results_HP <- Results_HP %>% 
  group_by(diffexpressed) %>%
  arrange(ifelse(diffexpressed == "DOWN", logFC, -logFC))%>%
  mutate(Top = ifelse(row_number() <= 10, Symbol,""))%>% 
  mutate(Top = ifelse(diffexpressed == "NO", "", Top)) 
row.names(Results_HP) <- Results_HP$ID

plot <- ggplot(data = Results_HP, aes(x = logFC, y = -log10(PValue), col = diffexpressed, label =Top)) +
  geom_point() +
  theme_test() +
  geom_text_repel(size=4, max.overlaps = 20, nudge_x = .1,
                  box.padding = 0.5,
                  point.padding = 0.01,
                  nudge_y = .3,
                  min.segment.length = 0,
                  segment.size= 0.5,
                  show.legend = FALSE) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "gray", "UP" = "red"),
                     breaks = c("DOWN", "UP"),
                     labels = c("DOWN", "UP")) +
  geom_vline(xintercept = -1, col = "gray30", linetype = "dashed") +
  geom_vline(xintercept = 1, col = "gray30", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.005337623), col = "gray30", linetype = "dashed") +
  theme(legend.title = element_blank(),axis.text=element_text(size=13), 
        axis.title=element_text(size=15),text =element_text(size=13))
png("Volcano_HP.png", units = "in", width=8, height=6, res=500)
plot
dev.off()

CPM_HP <- cpm(DGE_HP,normalized.lib.sizes = T, log = TRUE)

MinVals <- apply(CPM_HP, 1, min)
sum(MinVals == 0)
Exp <- CPM_HP[MinVals > 0, ]
Exp_log2 <- as.matrix(Exp)

png("boxplot_HP.png", units = "in", width=8, height=6, res=500)
boxplot(Exp_log2, ylab="log2 counts", main="Normalized Counts HP", las=1, cex.axis=0.8)
dev.off()

PCA <- prcomp(t(CPM_HP))
png("PCA_HP.png", units = "in", width=7, height=6, res=500)
autoplot(PCA,x=1, y=2, colour = "group", shape = "fraction",size=3,data = sample_table_HP) + 
  geom_hline(yintercept=0, col="gray") + geom_vline(xintercept=0,col="gray") + theme_minimal() + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) + 
  guides(colour=guide_legend(title="Group"),shape=guide_legend(title="Time"))
dev.off()

Cor_matrix <- cor(CPM_HP)
png("Cor_plot_HP.png", units = "in", width=7, height=6, res=500)
corrplot(Cor_matrix, method = 'color',is.corr = FALSE, col.lim = c(0.8, 1), col = COL1("Blues"))
dev.off()

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################


# pipeline comparisson 


desEQ_input <-  read.delim("DESeq.In.csv", sep = ",")
desEQ_Mono <-  read.delim("DESeq.Mon.csv", sep = ",")
desEQ_LP <-  read.delim("DESeq.LP.csv", sep = ",")
desEQ_HP <-  read.delim("DESeq.HP.csv", sep = ",")


colnames(desEq_HP) <- colnames(desEq_input)

Results_input$Gene <- rownames(Results_input)
Results_Mono$Gene <- rownames(Results_Mono)
Results_LP$Gene <- rownames(Results_LP)
Results_HP$Gene <- rownames(Results_HP)

length(intersect(desEq_Mono$Gene[desEq_Mono$Pvalue<0.05], 
                 Results_Mono$Gene[Results_Mono$PValue < 0.05]))

length(intersect(desEq_Mono$Gene, 
                 Results_Mono$Gene))



length(intersect(desEq_input$Gene[1:600], Results_input$Gene[1:600]))
table(desEq_Mono$Pvalue < 0.05)
table(Results_Mono$PValue < 0.05)


list_venn <- list(DownC1= desEQ_HP$X[desEQ_HP$P.value < 0.05 & desEQ_HP$log2.FC. < -1],
                  UpC1=  desEQ_HP$X[desEQ_HP$P.value < 0.05 & desEQ_HP$log2.FC. > 1],
                  DownC2= Results_HP$Gene[Results_LP$PValue < 0.05 & Results_HP$logFC < -1],
                  UpC2= Results_HP$Gene[Results_LP$PValue < 0.05 & Results_HP$logFC > 1])


ggvenn(list_venn)

intersect(desEQ_HP$X[desEQ_HP$P.value < 0.05 & desEQ_HP$log2.FC. < -1],
          Results_HP$Gene[Results_LP$PValue < 0.05 & Results_HP$logFC > 1])

rm(list = ls(pattern = "desEq"))
