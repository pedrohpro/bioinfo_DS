
#rm(list = ls())

# Library importation
library(tidyverse) #do it all
library(org.Hs.eg.db) #Homo sapiens OrgDb
library(biomaRt) # gene annotations
library(VennDiagram) #venn digrams
library(karyoploteR) #GWAS plots
library(ggrepel) #text repel

# for gene annotation
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # for gene annotation

# Data importation

#adip1 <- read.delim(file='data/vat.tsv')
#adip2 <- read.delim(file='data/sat.tsv')
#adipose <- rbind(adip1, adip2)
studies <- list()
studies$vat <- read.delim(file='data/vat.tsv')
studies$sat <- read.delim(file='data/sat.tsv')
studies$gut <- read.delim(file='data/gut.tsv')
studies$pleura <- read.delim(file='data/pleura.tsv')
studies$kidney <- read.delim(file='data/kidney.tsv')
studies$pbmc28 <- read.delim(file='data/pbmc28.tsv')

# Transforming celltype column into factor
studies <- map(studies, ~ mutate(., celltype = factor(celltype)))

# Adding/correcting columns
studies <- map(studies, function(df){
  
  # Features annotation
  features_annot <- getBM(
    attributes = c("external_gene_name", "chromosome_name", "start_position", "end_position"),
    filters = "external_gene_name",
    values = df$Genes,
    mart = ensembl
  )
  features_annot <- features_annot[features_annot$chromosome_name %in% c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"),]
  features_annot$chromosome_name <- factor(features_annot$chromosome_name, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"))
  features_annot$posit <- (features_annot$start_position + features_annot$end_position)/2
  
  # chromosome
  df$chr <- features_annot[match(df$Genes, features_annot$external_gene_name), "chromosome_name"]
  df$chr <- factor(df$chr)
  
  # position
  df$posit <- features_annot[match(df$Genes, features_annot$external_gene_name), "posit"]
  
  # Dropping features with no annotations
  df <- drop_na(df)
  
  # DEG
  df$DEG <- "NO"
  df$DEG[df$avg_log2FC > 0.585 & df$p_val_adj < 0.05] <- "UP"
  df$DEG[df$avg_log2FC < -0.585 & df$p_val_adj < 0.05] <- "DOWN"
  df$DEG <- factor(df$DEG)
  
  # color
  df$color[df$DEG == "NO"] <- "#cccccc"
  df$color[df$DEG == "UP"] <- "#FF8BA0" 
  df$color[df$DEG == "DOWN"] <- "#00AFBB"
  
  return(df)
})

# temporário pra mudar a cor sem rodar o time consuming biomart
# studies <- map(studies, function(df){
#   df$color[df$DEG == "NO"] <- "#cccccc"
#   df$color[df$DEG == "UP"] <- "#FF8BA0" 
#   df$color[df$DEG == "DOWN"] <- "#00AFBB"
#   
#   return(df)
# })


# List of lists by cell type
subsets <- list()
for (i in names(studies)) {
  study <- list()
  for (j in levels(studies[[i]]$celltype)) {
    study[[j]] <- studies[[i]][studies[[i]]$celltype == j,]
  }
  subsets[[i]] <- study
}

# Re-factoring cell type columns after subsetting and adding gene highlight column
subsets <- map(subsets, function(i){
  i <- map(i, function(j){
    # Refactoring
    j$celltype <- factor(j$celltype)
    
    # Highlight genes
    aux_df <- j[j$p_val_adj < 0.05 & abs(j$avg_log2FC) < 3,]
    aux_df <- aux_df[order(-abs(aux_df$avg_log2FC)),]
    mark <- aux_df[1:ifelse(nrow(aux_df) < 40, nrow(aux_df), 40),]$Genes
    j$highlight <- NA
    j$highlight[j$Genes %in% mark] <- j[j$Genes %in% mark,]$Genes
    
    return(j)
  })
  return(i)
})

# Plotting manhattan and karyoplot
for (stu_name in names(subsets)) {
  for (cell_df in subsets[[stu_name]]) {
    ########################
    #####  Manhattan  ######
    ########################
    annot_df <- cell_df %>%
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
    
    # PLOT
    manhat <- ggplot(cell_df, aes(x = chr, y = avg_log2FC, color = DEG, label = highlight))+
      geom_jitter(size = 0.6)+
      labs(
        title = paste0("DEGs across chromosomes for ", levels(cell_df$celltype), " from ", stu_name),
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
      geom_text_repel(
        max.overlaps = Inf, 
        color = cell_df$color, ##### 
        size = 2
      )+ # To show all labels
      annotate(
        'rect',
        xmin = annot_df$x_axis - 0.4,
        xmax = annot_df$x_axis + 0.4,
        ymin = annot_df$y_axis + 0.3,
        ymax = annot_df$y_axis - 0.3,
        alpha = 0.6,
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
    
    ggsave(paste0("results/manhat_", levels(cell_df$celltype), "_", stu_name, ".png"), plot = manhat, height = 5, width = 10)
  
    ########################
    #####  Karyoplot  ######
    ########################
    
  }
}










############################
#####  Venn diagrams  ######
############################

# Getting DEGs by study
degs <- list()
for (i in names(subsets)) { # iterate through studies names
  for (j in names(subsets[[i]])) { # iterate through cell types
    degs[[i]][[j]] <- subsets[[i]][[j]][subsets[[i]][[j]]$DEG != "NO", ]$Genes
  }
}


















########################
#####  Karyoplot  ######
########################


# Preparing the df
local_res_df <- i[i$chr == "X",]
local_res_df$chr <- paste0("chr", local_res_df$chr)
View(local_res_df)

# Selecting genes to highlight
#aux1_df <- local_res_df[local_res_df$alterat != "NO",]
#aux1_df <- aux1_df[order(-abs(aux1_df$log2FoldChange)),]
#aux1_df <- aux1_df[1:5,] # top 5 highest signif log2fc
#aux2_df <- local_res_df[rownames(local_res_df) %in% intGenes,] # genes of interest from chr21
#local_highlight_df <- rbind(aux1_df, aux2_df)

# Initializing save plot
png(filename = paste0("results/karyoplot_", levels(i$celltype), ".png"), width = 1500, height = 500)

# Plotting
kp <- plotKaryotype(genome="hg38", chromosomes="chrX", main = paste0("Gene alterations in the X chromosome for ", levels(i$celltype)), cex = 2)

kpAddCytobandLabels(kp, 
                    force.all = TRUE, #all name sections of the chromosome
                    srt = 90, 
                    cex = 1, 
                    col = c(rep("black", times = 7), "white", "black", "white", rep("black", times = 12), "white", "black", "white",  rep("black", times = 7), "white", rep("black", times = 5), "white", "black") #to contrast the backgorund better
)

# kpPoints(kp, 
#          data = local_res_df, 
#          chr = "chrX", #local_res_df$chr, 
#          x = local_res_df$posit, 
#          y = local_res_df$avg_log2FC, 
#          ymax = 7, #max(local_res_df$log2FoldChange), 
#          ymin = -7, #min(local_res_df$log2FoldChange),
#          #col = as.character(local_res_df$color),
#          cex = 1.8
# )

# kpAxis(kp, 
#   ymax = 3, #max(local_res_df$log2FoldChange), 
#   ymin = -3, #min(local_res_df$log2FoldChange),
#   tick.pos = NULL, #c(0, max(local_res_df$log2FoldChange), min(local_res_df$log2FoldChange), max(local_res_df$log2FoldChange)/2, min(local_res_df$log2FoldChange)/2),
#   numticks = 7, 
#   cex = 0.9
# )
# 
# kpAbline(kp, h = 0.5) # this might to be adjusted with kpAxis()
# 
# kpAddLabels(kp, 
#   labels = "log2 Fold Change", 
#   srt = 90, 
#   pos = 1, 
#   label.margin = 0.04,
#   cex = 1.5
#   #ymax = max(local_res_df$log2FoldChange), 
#   #ymin = min(local_res_df$log2FoldChange)
# )

#kpPlotMarkers(kp, 
#  chr = "chr21",
#  x = local_highlight_df$position, 
#  y = rescale(local_highlight_df$avg_log2FC, from = c(-3, 3), to = c(-4, 1)),
#  labels = rownames(local_highlight_df), 
#  cex = 1.2, 
#  r0 = 0.8,
#  line.color = NULL,
#  text.orientation = "vertical",
#  adjust.label.position = TRUE
#)

# Closing save plot
dev.off()