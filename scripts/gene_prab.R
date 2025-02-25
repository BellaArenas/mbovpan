# Author: Noah Legall
# Purpose: Output a presence absence matrix of M. bovis virulence genes 
library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)
library(ggdendro)
library(ggtree)
library(ggnewscale)

# load in the general metadata
args = commandArgs(trailingOnly=TRUE)
gene_pres_abs <- read.csv(args[1], header = TRUE, stringsAsFactors = FALSE, row.names = "Gene")
print("the prab was read in successfully")
# load in the gene presence absence data, keep only accessory
# we should already have access to this in our directory
accessory_pa <- gene_pres_abs %>% select(3:(ncol(gene_pres_abs)))


accessory_pa[!(accessory_pa=="")] <- 1
accessory_pa[accessory_pa==""] <- 0

print("The dimensions of accessory_pa")
print(ncol(accessory_pa))
print(nrow(accessory_pa))

num_col <- ncol(accessory_pa)

task <- function(x){sum(as.numeric(as.character(x)))}

accessory_pa$pr <- apply(accessory_pa,1,task)
accessory_pa$perc_pr <- accessory_pa$pr/num_col
accessory_pa <- accessory_pa %>% filter(perc_pr >= 0.15 & perc_pr <= 0.99) %>% select(-c("pr","perc_pr"))

print("percentages are calculated")
head(accessory_pa)
nrow(accessory_pa)

accessory_pa$gene_id = rownames(accessory_pa)


accessory_pa_long <- accessory_pa %>% gather(sample,prab,-gene_id)
accessory_matrix <- as.matrix(accessory_pa %>% select(-gene_id))

# perform clustering 
accessory_transpose <- t(accessory_matrix)
rownames(accessory_transpose) <- gsub("_trimmed_R1.scaffold.annot","",rownames(accessory_transpose))


accessory_matrix <- as.matrix(accessory_pa %>% select(-gene_id))
rownames(accessory_matrix) <- accessory_pa$gene_id
accessory_dendro <- as.dendrogram(hclust(d = dist(accessory_transpose, method = "binary"), method = "ward.D"))

ad_gg <- ggtree(accessory_dendro)
ad_gg[["data"]]$label <- gsub(".annot","",ad_gg$data$label)

# check if the dendrogram was made accordingly based on the data present
#print(ad_gg[["data"]])

pdf("gene_prab_figures.pdf")

 isolate_dat <- read.csv(args[2], stringsAsFactors = FALSE, check.names = FALSE)
 ad_gg <- ad_gg %<+% isolate_dat +
        ggtree::vexpand(.1, -1)

 
#function to make the dendrograms programatically
    mbov_tree <- function(mytree,metadata){
      

      #checking how the data looks after incorporationg metadata

      

      row_id <- subset(mytree[["data"]], isTip == TRUE)$label
      mytree_onlytip <- as.data.frame(subset(mytree[["data"]], isTip == TRUE)[,metadata])
      rownames(mytree_onlytip) <- row_id



      t1 <- gheatmap(mytree, mytree_onlytip, width = 0.3, colnames = FALSE) +
        scale_fill_manual(values = hcl.colors(length(unique(mytree_onlytip[,metadata])),palette = "Zissou 1"), name = metadata)
      
      #plot(t1)
      t1_scaled <- t1 + new_scale_fill()
      t2 <- gheatmap(t1_scaled,accessory_transpose, offset = 3, colnames_angle=45,hjust = 1,colnames_offset_y = -2.5,font.size = 3) +  
        scale_fill_manual(values = c("gray75","darkblue"), name = "Presence/Absence")
      Sys.sleep(1) #gives program the time to make the figure
      print(t2[["data"]])
      return(t2)
      
      }




 
  # this works out 
  #head(isolate_dat)
  print(colnames(isolate_dat))
  for(i in 1:length(colnames(isolate_dat))){
    print(colnames(isolate_dat)[i])
    if(colnames(isolate_dat)[i] == "Name"){ 
      print("skipping the Name column")
      next
    }

    else if(length(unique(isolate_dat[,i])) < 2){
      print("not unique metadata")
      next
    }
    else{
      
      print(mbov_tree(ad_gg,colnames(isolate_dat)[i]))
    
  
} 
  }

dev.off()

### Finish Code and update mbovpan 