library(tidyverse)
library(readxl)
library(biomaRt)
library(DESeq2)
library(geneEnrichment)

modelgenes <- read_csv("../1_CriGr_model_ref_Data/genes_mapped.csv")
rnaseq <- read_csv("data/reads_concatenated_rmvD5.csv") %>% 
  filter(!if_all(matches("Reads"), ~ .x < 10))

### SET UP BIOMART
# Must be connected to the internet for this!
ensembl_cho <- useEnsembl(biomart = "genes",
                          dataset = "cgchok1gshd_gene_ensembl")
ensembl_mm <- useEnsembl(biomart = "genes",
                         dataset = "mmusculus_gene_ensembl")
# I've commented out some common identifers/ attributes. 
# Change according to your needs! :)
bm_dat_cho <- getBM(attributes = c(#"uniprotswissprot",
  #  "uniprot_gn_symbol",
  "ensembl_gene_id",
  #  "ensembl_transcript_id",
  "entrezgene_id"),
  filters = "ensembl_gene_id",
  values = rnaseq$gene, 
  mart = ensembl_cho) %>% 
  dplyr::select(GeneID= entrezgene_id, Ensembl = ensembl_gene_id) 

# Wrangle data into DESeq2 format.
rnaseq_forDeSeq <- rnaseq %>% 
  tibble::column_to_rownames("gene") %>% 
  dplyr::select(matches("^B"))

colnames(rnaseq_forDeSeq) <- gsub("_ReadsPerGene\\.out\\.tab", "", colnames(rnaseq_forDeSeq))

# All the reactors were "high glutamine" condition.
# We have 2 technical reps. per reactor and 3 biological reps.
# All these DESeq commands were adapted from the documentation!
sampleData <- data.frame(
  sampleName = colnames(rnaseq_forDeSeq), 
  time = gsub("B1[234]_D", "", 
              gsub("_[12]$", "", colnames(rnaseq_forDeSeq))
  ),
  row.names = colnames(rnaseq_forDeSeq)
) %>% 
  arrange(as.integer(time)) #Rearrange in time order

#Rearrange in time order
rnaseq_forDeSeq <- rnaseq_forDeSeq[,sampleData$sampleName]

dds <- DESeqDataSetFromMatrix(countData = rnaseq_forDeSeq,
                              colData = sampleData,
                              design =  ~ time)
dds <- DESeq(dds, fitType = "local")
rld <- rlog(dds)
pca <-  plotPCA(rld, intgroup = "time",ntop= length(rld), returnData = T)
pca$group <- factor(as.numeric(as.character(pca$group)))
pca$rep <- ifelse(grepl("B12", pca$name), 1,
                  ifelse(grepl("B13", pca$name), 2,3))

palette = colorRampPalette(c("red", "gold", "darkgreen", "darkblue","violet"))(15)
ggplot(pca, aes(x = PC1, y = PC2, color = group, shape = as.factor(rep))) +
  geom_point(size = 4) +
  xlab("PC1 61.9%") +
  ylab("PC2 20.0%") +
  scale_color_manual(name = "Day",
                   #  values = colorRampPalette(c("red", "gold", "darkgreen", "darkblue","violet"))(7)
                     values = c(
                       palette[5],
                       palette[7],
                       palette[8],
                       palette[9],
                       palette[12],
                       palette[13],
                       palette[15]
                     )
  )
plotDispEsts(dds)

# Set up contrasts for stats tests
# Concurrent days were contrasted 
# (i.e. day 4 - day 6, day 6 - day 7, etc.)
timeComp1 <- unique(sampleData$time)[1:length(unique(sampleData$time))-1]
timeComp2 <- unique(sampleData$time)[2:length(unique(sampleData$time))]
# Perform tests. All parameters default.
res <- lapply(1:length(timeComp1), function(i){
  #res <- 
  results(dds, contrast=c("time", timeComp2[i], timeComp1[i])) %>% 
    as.data.frame() %>% 
    rownames_to_column("gene") %>% 
    dplyr::select(gene, log2FoldChange, padj) %>% 
    mutate(contrast  = paste0(timeComp1[i], "v", timeComp2[i]),
           time = timeComp2[i])
  #return(res[res$padj < 0.05 & !is.na(res$padj),])
})

# Significantly changing genes = FDR < 0.05 for at least one time point contrast
sig_genes <- res %>% 
  bind_rows() %>% 
  group_by(gene) %>% 
  filter(padj < 0.05 & !is.na(padj))
length(unique(sig_genes$gene))
table(modelgenes$Ensembl %in% sig_genes$gene) # n model genes differentially expressed over time

# Extract normalized counts from DeSeq2 object.
ntd <- counts(dds, normalized = TRUE)

# Tidy data for export
sig_genes_forExport <- sig_genes %>% 
  dplyr::select(-time) %>% 
  pivot_longer(c(padj, log2FoldChange)) %>% 
  mutate(name = paste0(name, "_", contrast)) %>% 
  dplyr::select(-contrast) %>% 
  pivot_wider()
ntd %>% 
  as.data.frame() %>% 
  rownames_to_column("Ensembl") %>% 
  left_join(sig_genes_forExport, by = c("Ensembl" = "gene")) %>%  
  left_join(bm_dat_cho, by = "Ensembl", multiple = "all") %>% 
  mutate(isSig = ifelse(Ensembl %in% sig_genes$gene, "+", "")) %>% 
  write_csv("results/data/DeSeq2_counts_wDE_genes_251124_rmvD5.csv")

# Calculate Z-score for heatmap visualisation
z_long <-  ntd %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  filter(gene %in% sig_genes$gene) %>% 
  pivot_longer(-gene,
               names_to = "day") %>% 
  mutate(day = as.integer(gsub("_", "", substr(day, 6,7)))) %>% 
  group_by(gene,day) %>% 
  summarise(m = mean(value, na.rm = T), .groups = "keep") %>% 
  ungroup() %>% 
  group_by(gene) %>% 
  mutate(z = (m - mean(m)) / sd(m))

# Data are Z normalised and ordered along the y axis by the time point of peak expression.
z_long %>% 
  group_by(gene) %>% 
  filter(sum(is.na(z)) < n()) %>% 
  mutate(max_tp = factor(day[z == max(z)],
                         levels = c(4,6,7,8,11,12, 14),
                         ordered = T)) %>% 
  ungroup() %>% 
  arrange(max_tp) %>% 
  mutate(gene =factor(gene, levels = unique(gene))) -> z_long2

z_long2 %>% 
  ggplot(aes(x = as.factor(day),
             y = gene,
             fill = z)) +
  geom_tile() +
  scale_fill_gradient2(high = "red", mid = "white", low = "blue") +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank()) +
  xlab("Day")  +
  ggtitle("DE genes") -> z_g
z_g

z_g_model <- z_g$data %>% 
  filter(gene %in% modelgenes$Ensembl) %>% 
  ggplot(aes(x = as.factor(day),
             y = gene,
             fill = z)) +
  geom_tile() +
  scale_fill_gradient2(high = "red", mid = "white", low = "blue") +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank())+
  xlab("Day") +
  ggtitle("",
          "filtered for genes in GEM (1158 / 2172)") 
# To check # genes in GEM (multiple Ensembl IDs for some NCBI Gene IDs): 
# z_g_model$data %>%  left_join(modelgenes, by = c("gene" = "Ensembl")) %>% dplyr::select(GeneID) %>%  distinct()
z_g_model

# Visualise n. genes in each sample
ntd %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  filter(gene %in% sig_genes$gene) %>% 
  pivot_longer(-gene,
               # names_to = "day"
  ) %>% 
  mutate(day = as.integer(gsub("_", "", substr(name, 6,7)))) %>% 
  filter(value > 10) %>% 
  group_by(day, name) %>% 
  summarise(n = n()) %>%
  group_by(day ) %>% 
  summarise(mean_n = mean(n), sd_n = sd(n))


# To calculate enrichment of pathways in each "cluster"
mm_ensembl <- getBM(attributes = c(
  #"uniprot_gn_symbol",
  "mmusculus_homolog_ensembl_gene",
  "ensembl_gene_id"),
  filters = "ensembl_gene_id",
  values = z_long$gene, 
  mart = ensembl_cho)  
mm_symbols <- getBM(attributes = c(
  #"uniprot_gn_symbol",
  "mgi_symbol",
  "ensembl_gene_id"),
  filters = "ensembl_gene_id",
  values = mm_ensembl$mmusculus_homolog_ensembl_gene, 
  mart = ensembl_mm)  

enrich <-  left_join(mm_ensembl, mm_symbols, by = c("mmusculus_homolog_ensembl_gene" = "ensembl_gene_id"), multiple= "all") %>% 
  left_join(z_long2, by = c("ensembl_gene_id" = "gene")) %>%
  filter(if_all(c(max_tp), ~ !is.na(.x))) %>% 
  #  mutate(timePhase = max_tp) %>% 
  mutate(timePhase = ifelse(max_tp ==4, "4",
                            ifelse(max_tp %in% c(6,7,8), "6,7,8", "11,12,14"))) %>% 
  group_by(timePhase) %>% 
  group_split() %>% 
  map(~ calculateEnrichment(unique(.x$mgi_symbol),
                            "KEGG_2021_Mouse",
                            #"GO_Biological_Process_2023", 
                            visualise = F, simplify = T) %>% 
        mutate(timePhase = unique(.x$timePhase)))

enrich_g <- enrich %>% 
  map(~ .x %>% 
        arrange(Adjusted.P.value) %>% 
        dplyr::slice(1:10)
  ) %>%  
  bind_rows() %>% 
  group_by(Term) %>% 
  mutate(n = n()) %>% 
  group_by(n) %>% 
  arrange(timePhase) %>% 
  ungroup() %>% 
  mutate(n = as.integer(gsub("\\/.*", "", Overlap)),
         order = 1:n()) %>% 
  ggplot(aes(x = as.factor(timePhase), y = reorder(Term, order))) +
  geom_point(aes(color = as.factor(timePhase), size = n), alpha = 0.7) +
  xlab("Cluster") +
  ylab("") +
  guides(color ="none") +
  scale_color_manual(values = c(
    "4"="red",
    "6,7,8"="darkgreen", 
    "11,12,14" = "darkblue"
  )) +
  scale_x_discrete(limits = c("4", "6,7,8", "11,12,14"),
                   breaks = c("4", "6,7,8", "11,12,14")) +
  # scale_x_discrete(breaks = factor(c(4,6,7,8,11,12,14)),
  #                  limits= factor(c(4,6,7,8,11,12,14)))
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()

enrich_g


enrich_flt <- enrich %>% 
  map(~ .x %>% 
        arrange(Adjusted.P.value) %>% 
        dplyr::slice(1:10)
  ) %>%  
  bind_rows() %>% 
  dplyr::select(Term, timePhase, Genes) %>% 
  group_by(Term) %>% 
  mutate(isUnique = ifelse(n() == 1, "+","") ) %>% 
  write_csv("results/data/timePhase_enrichment_251124_rmvD5.csv")
