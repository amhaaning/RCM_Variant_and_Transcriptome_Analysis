library(ggplot2)
library(networkD3)
library(dplyr)
library(viridisLite)
library(scales)

setwd(".")

#Make a Sankey diagram - connects genes with their associated GO terms####

#Read in data for Sankey plot (top 50 genes based on p-value)
links <- read.delim("./topGO_Enrichment_Results/topGO_edgeR_results_for_Sankey_plot.txt")

links <- links[,c(3,1,4,5)]

links$value <- 1

colnames(links)[1:2] <- c("source","target")

nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

nodes$group <- c("0","1","2","3","4","5","6","7",paste(rep("gene",31),sep=","))

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

links$group <- factor(links$IDsource)

#Set colors for links and nodes
my_color <- 'd3.scaleOrdinal() .domain(["0","1","2","3","4","5","6","7" ,"gene"]) .range(["#33CCFF","#FFCCFF","#FF9999","#66FFFF","#FFFF99","#9999FF","#FF66CC","#66FF66","#5F5F5F"])'

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource",Target = "IDtarget", Value = "value", NodeID = "name", colourScale=my_color, NodeGroup="group", LinkGroup="group", nodeWidth = 20,fontSize = 20,sinksRight=FALSE,margin = list(top=4, right=2, bottom=4, left=2))

saveNetwork(p, "sankey_diagram_GO_Analysis.html", selfcontained = TRUE)

#Visualize results of topGO enrichment####

#Read in topGO results
topGO_MF <- read.delim("./topGO_Enrichment_Results/All_DE_Genes/RCM_All_DE_Genes_MF_topGO_enrichment.txt")
topGO_CC <- read.delim("./topGO_Enrichment_Results/All_DE_Genes/RCM_All_DE_Genes_CC_topGO_enrichment.txt")

#Add GO category
topGO_MF$Category <- "Molecular\nFunction"
topGO_CC$Category <- "Cellular\nCompartment"

#Combine all signficantly associated GO terms
topGO <- rbind(topGO_MF,topGO_CC)

#Sort by p-value and keep the top 25 associated terms
topGO <- topGO[order(topGO$weightFisher_adj),]
topGO <- topGO[1:25,]

#Make a column with log10 pvalue
topGO$log10pval <- -1*(log10(topGO$weightFisher))

#Make a column with fold-difference (observed/expected)
topGO$fold_diff <- topGO$Significant/topGO$Expected

#Make a column with GO term and description
topGO$GO_Term <- paste0(topGO$GO.ID,": ",topGO$Term)

#Order GO terms by p-value
topGO$GO_Term <- factor(topGO$GO_Term,levels=topGO$GO_Term[order(topGO$log10pval)])

mypal <- colorRampPalette(c("navy", "white", "red"))(20)

#Make overrepresentation plot
pdf("Overrepresentation_plot_top_25_GO_terms.pdf",height=5,width=8)
ggplot(topGO,aes(x=log10pval,y=GO_Term,fill=fold_diff)) + geom_bar(stat = "identity") + theme_bw() + facet_grid(rows="Category",space = "free_y",scales="free_y") + theme(axis.ticks=element_blank(),axis.title.y=element_blank(),axis.text=element_text(color="black",size = 11),axis.title.x=element_text(color="black",size=12),legend.title=element_text(color="black",size=11),legend.text = element_text(color="black",size=9),strip.text=element_text(size=11,color="black"),legend.title.align = 0.5) + xlab("-log10(weighted fisher p-value)") + scale_fill_gradientn(colors=c("#404788FF","#1F968BFF","#73D055FF"),values=rescale(c(1.09,1.5,9.804)),limits=c(1.09,9.804),breaks=c(1.25,2,3,4,9)) + scale_x_continuous(limits=c(0,12),expand=c(0,0)) + labs(fill="fold-difference\n(observed/expected)")
dev.off()
