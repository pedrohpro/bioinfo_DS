#rm(list = ls())

# Library importation
library(tidyverse) #do it all
library(org.Hs.eg.db) #Homo sapiens OrgDb
library(biomaRt) # gene annotations
library(karyoploteR) #GWAS plots
library(ggrepel) #text repel

# for gene annotation
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # for gene annotation

# Data importation
df <- read.delim(file='GSE190125.top.table.tsv')

# Features annotation
features_annot <- getBM(
  attributes = c("external_gene_name", "chromosome_name", "start_position", "end_position"),
  filters = "external_gene_name",
  values = df$Symbol,
  mart = ensembl
)

# getting only chromosome genes
features_annot <- features_annot[features_annot$chromosome_name %in% c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"),]

# re-leveling factor to order chromosomes
features_annot$chromosome_name <- factor(features_annot$chromosome_name, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"))

# Calculating avg gene position
features_annot$posit <- (features_annot$start_position + features_annot$end_position)/2
  
# adding chromosome annotation
df$chr <- features_annot[match(df$Symbol, features_annot$external_gene_name), "chromosome_name"]
df$chr <- factor(df$chr)
  
# adding position annotations
df$start_position <- features_annot[match(df$Symbol, features_annot$external_gene_name), "start_position"]
df$end_position <- features_annot[match(df$Symbol, features_annot$external_gene_name), "end_position"]
df$avg_position <- features_annot[match(df$Symbol, features_annot$external_gene_name), "posit"]
  
# Dropping features with no annotations
df <- drop_na(df)
  
# DEG
df$DEG <- "NO"
df$DEG[df$log2FoldChange > 1 & df$padj < 0.05] <- "UP"
df$DEG[df$log2FoldChange < -1 & df$padj < 0.05] <- "DOWN"
df$DEG <- factor(df$DEG)
  
# color
#df$color[df$DEG == "NO"] <- "#cccccc"
#df$color[df$DEG == "UP"] <- "#bb0c00" 
#df$color[df$DEG == "DOWN"] <- "#00AFBB"
  
# Highlight genes
# aux_df <- df[df$padj < 0.05,] #& abs(df$log2FoldChange) < 3
# aux_df <- aux_df[order(-abs(aux_df$log2FoldChange)),]
# mark <- aux_df[1:ifelse(nrow(aux_df) < 40, nrow(aux_df), 40),]$Symbol
# df$highlight <- NA
#df$highlight[df$Symbol %in% mark] <- df[df$Symbol %in% mark,]$Symbol

# list of dfs by chromosome
chr_list <- list()
for (i in levels(df$chr)) {
  aux_df <- df[df$chr == i,]
  aux_df$chr <- factor(aux_df$chr)
  aux_df <- aux_df[order(aux_df$start_position),]
  chr_list[[i]] <- aux_df
}

# clustering loop
chr_list <- map(chr_list, function(aux){
  cluster <- 1
  aux$cluster <- NA
  aux$cluster[1] <- cluster
  for (i in 2:nrow(aux)) {
    if (aux$end_position[i-1] - aux$start_position[i] < 1E5){ # how many clusters do you want??? | aux$end_position[i-1] > aux$start_position[i]
      aux$cluster[i] <- cluster
    }else{
      cluster <- cluster + 1
      aux$cluster[i] <- cluster
    }
  }
  aux$cluster <- factor(aux$cluster)
  return(aux)
})

# list of chr separated by cluster
cluster_list <- list()
for (i in names(chr_list)) {
  aux <- list()
  for (j in levels(chr_list[[i]]$cluster)){
    aux[[j]] <- chr_list[[i]][chr_list[[i]]$cluster == j,]
    aux[[j]]$cluster <- factor(aux[[j]]$cluster)
  }
  cluster_list[[i]] <- aux
}

# Cluster df by chromosome
chr_clusters <- list()
for (i in names(chr_list)) {
  chr_clusters[[i]] <- chr_list[[i]] %>% 
    group_by(cluster) %>% 
    summarise(num_notdeg = sum(DEG == "NO"), num_degs = sum(DEG != "NO")) %>% 
    ungroup()
}

# include other squares of the contingency matrix
chr_clusters <- map(chr_clusters, function(i){
  i$num_outdegs <- nrow(df[df$DEG != "NO",]) - i$num_degs
  i$num_outnotdeg <- nrow(df[df$DEG == "NO",]) - i$num_notdeg
  return(as.data.frame(i))
})

# include cluster location
for(chr in names(chr_clusters)){
  chr_clusters[[chr]]$start_position <- NA
  chr_clusters[[chr]]$end_position <- NA
  chr_clusters[[chr]]$chr <- paste0("chr", chr)
  for (num_cluster in names(cluster_list[[chr]])) {
    chr_clusters[[chr]]$start_position[chr_clusters[[chr]]$cluster == num_cluster] <- min(cluster_list[[chr]][[num_cluster]]$start_position)
    chr_clusters[[chr]]$end_position[chr_clusters[[chr]]$cluster == num_cluster] <- max(cluster_list[[chr]][[num_cluster]]$end_position)
  }
}

# chi-square test
calculate_chi2 <- function(row_data) {
  observed <- matrix(row_data, nrow = 1)
  expected <- rep(sum(row_data) / length(row_data), length(row_data))
  result <- chisq.test(observed, p = expected)
  p_value <- result$p.value
  return(p_value)
}
# Apply the function to each row of the data frame
chr_clusters <- map(chr_clusters, function(i){
  i$pvalue <- NA
  i$pvalue <- apply(i[,c("num_degs","num_outdegs","num_notdeg","num_outnotdeg")], 1, calculate_chi2) #NOTE: order of columns!!!
  return(i)
})
# ???
chr_clusters <- map(chr_clusters, function(i){
  i$pvalue <- NA
  for (clu in i$cluster) {
    chisq <- chisq.test(matrix(i[clu,c("num_notdeg","num_degs","num_outnotdeg","num_outdegs"), drop = TRUE], nrow = 2, byrow = TRUE))
    i$pvalue[i$cluster == clu] <- chisq[["p.value"]]
  }
  return(i)
})

# calculate padj by "BH" and enrichment score
chr_clusters <- map(chr_clusters, function(i){
  i$padj <- p.adjust(i$pvalue, method = "BH")
  i$score <- -log10(i$padj)
  return(i)
})

# karyoplot with density/histogram plot + scatter por cima pra mostrar as degs?
############# JUST A TEST!!!!!
chr_clusters <- map(chr_clusters, function(i){
  i$score <- sample(0:6, nrow(i), replace = TRUE)
  return(i)
})

# combining all clusters from all chromosomes
all_clusters <- as.data.frame(do.call(rbind, chr_clusters))

# Initializing save plot
png(filename = "results/karyoplot.png", width = 5000, height = 6000)

# Plotting
kp <- plotKaryotype(genome="hg38", main = "Genome wide DEG clusters density", cex = 7)

kpAddCytobandLabels(
  kp, 
  force.all = TRUE, #all name sections of the chromosome
  srt = 90, 
  cex = 1, 
  #col = c(rep("black", times = 6), "white", rep("black", times = 7)) #to contrast the background better
)

# kpPoints(kp, 
#          data = local_res_df, 
#          chr = local_res_df$chr, 
#          x = local_res_df$position, 
#          y = local_res_df$log2FoldChange, 
#          ymax = 3, #max(local_res_df$log2FoldChange), 
#          ymin = -3, #min(local_res_df$log2FoldChange),
#          col = as.character(local_res_df$color),
#          cex = 1.8
# )
# 

kpAxis(
  kp, 
  ymax = max(all_clusters$score),
  ymin = 0, #min(local_res_df$log2FoldChange),
  tick.pos = NULL, #c(0, max(local_res_df$log2FoldChange), min(local_res_df$log2FoldChange), max(local_res_df$log2FoldChange)/2, min(local_res_df$log2FoldChange)/2),
  numticks = 7, 
  cex = 0.9
)
# 
# kpAbline(kp, h = 0.5) # this might to be adjusted with kpAxis()
# 
kpAddLabels(kp, 
  labels = "enrichment score", 
  srt = 90, 
  pos = 4, 
  label.margin = 0.04,
  cex = 1.5
  #ymax = max(local_res_df$log2FoldChange), 
  #ymin = min(local_res_df$log2FoldChange)
)
# 
# kpPlotMarkers(kp, 
#               chr = "chr21",
#               x = local_highlight_df$position, 
#               y = rescale(local_highlight_df$log2FoldChange, from = c(-3, 3), to = c(-4, 1)),
#               labels = rownames(local_highlight_df), 
#               cex = 1.2, 
#               r0 = 0.8,
#               line.color = NULL,
#               text.orientation = "vertical",
#               adjust.label.position = TRUE
# )

# kpSegments(kp, 
#   chr = "chr21", 
#   x0 = local_highlight_df$position,
#   x1 = local_highlight_df$position,
#   y0 = local_highlight_df$log2FoldChange,
#   y1 = max(local_res_df$log2FoldChange),
#   ymax = max(local_res_df$log2FoldChange),
#   ymin = min(local_res_df$log2FoldChange),
#   r0 = 0.8
# )

# Closing save plot
dev.off()



















# Plotting manhattan and karyoplot
for (study_df in studies) {
  ########################
  #####  Manhattan  ######
  ########################
  annot_df <- study_df %>%
    group_by(chr) %>%
    summarise(
      up_degs = sum(DEG == "UP"), 
      down_degs = sum(DEG == "DOWN")
      #up_pct = (sum(DEG == "UP"))/(sum(DEG == "UP") + sum(DEG == "DOWN")),
      #down_pct = (sum(DEG == "DOWN"))/(sum(DEG == "UP") + sum(DEG == "DOWN"))
    ) %>%
    pivot_longer(cols = c("up_degs", "down_degs"), names_to = "cat", values_to = "num_degs")
  annot_df$num_degs <- as.character(annot_df$num_degs)
  
  # x column
  annot_df$x_axis <- NA
  for (i in 1:nlevels(annot_df$chr)) {
    annot_df$x_axis[annot_df$chr == levels(annot_df$chr)[i]] <- i
  }
  # color column
  annot_df$color <- NA
  annot_df$color[annot_df$cat == "up_degs"] <- "#bb0c00" #"#ffcccb"
  annot_df$color[annot_df$cat == "down_degs"] <- "#0044aa" #"#add8e6"
  # y column
  annot_df$y_axis <- NA
  annot_df$y_axis[annot_df$cat == "up_degs"] <- 3
  annot_df$y_axis[annot_df$cat == "down_degs"] <- -3
  
  manhat <- ggplot(study_df, aes(x = chr, y = log2FoldChange, color = DEG, label = highlight))+
    geom_jitter(size = 0.6)+
    labs(
      title = paste0("DEGs across chromosomes for PBMC in ", "Down Syndrome"),
      x = "Chromosome",
      y = "Log2 Fold Change",
      color = NULL
    )+ 
    # zoom
    coord_cartesian(ylim = c(-3, 3)) + # this puts genes outside these limits outside the plot
    scale_color_manual(values = c("DOWN" = "#0044aa77", "NO" = "#bbbbbb77", "UP" = "#bb0c0077")) +
    theme(   
      panel.background = element_rect(fill = "white"),
      panel.grid = element_line(color = "gray", linetype = "dotted"),
      axis.text.x = element_text(size = 12, face = "bold")  # Modify font size and make it bold
    )+
    geom_text_repel(max.overlaps = Inf, color = "black", size = 2)+ # To show all labels
    annotate(
      'rect',
      xmin = annot_df$x_axis - 0.4,
      xmax = annot_df$x_axis + 0.4,
      ymin = annot_df$y_axis + 0.3,
      ymax = annot_df$y_axis - 0.3,
      alpha = 0,
      fill = "white", #'#0044aa',
      col = "transparent" # Border color parameter
    ) +
    annotate(
      geom = "text",
      x = annot_df$x_axis, 
      y = annot_df$y_axis, 
      label = annot_df$num_degs,
      color = annot_df$color,
      size = 2.5
    )
  
  ggsave(paste0("manhatt_", ".png"), plot = manhat, height = 5, width = 10)
  
  ########################
  #####  Karyoplot  ######
  ########################
  
}
