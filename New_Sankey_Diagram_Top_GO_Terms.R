library(dplyr)
library(networkD3)

#Set working directory
setwd(".")

#Read in top 10 GO terms (fold-enrichment > 1.5 and then ranked by p-value)
topGO_Terms <- read.delim("Top10_GO_Terms_pval_and_fold_enrichment_1.5.txt",header=F)

#Read in all DE genes associated with significantly overrepresented GO terms
CC_GO <- read.delim("All_DE_Genes/RCM_All_DE_Genes_CC_associated_DE_Genes.txt")
MF_GO <- read.delim("All_DE_Genes/RCM_All_DE_Genes_MF_associated_DE_Genes.txt")

#Keep genes associated with top go terms
CC_GO <- CC_GO[CC_GO$Term %in% topGO_Terms$V1,]
MF_GO <- MF_GO[MF_GO$Term %in% topGO_Terms$V1,]

#Read in DE gene results from EdgeR
DE_Genes <- read.delim("../EdgeR_results_qlf_test_all_genes.txt")

#Keep gene symbol and p-values
DE_Genes <- DE_Genes[,c("GeneSymbol","PValue")]
colnames(DE_Genes) <- c("Associated_DE_Gene","EdgeR_pval")

#Merge genes associated with top GO terms with p-values from EdgeR analysis
CC_GO <- merge(CC_GO,DE_Genes,by="Associated_DE_Gene")
MF_GO <- merge(MF_GO,DE_Genes,by="Associated_DE_Gene")

#Keep top 5 genes by EdgeR p-value that are associated with the top GO terms
CC_GO_Sub <- CC_GO %>%
  arrange(EdgeR_pval) %>%
  group_by(GO.ID) %>%
  slice(1:6)

MF_GO_Sub <- MF_GO %>%
  arrange(EdgeR_pval) %>%
  group_by(GO.ID) %>%
  slice(1:6)

#Reformat results for Sankey diagram
CC_GO_Sub$GO_Term <- paste0(CC_GO_Sub$GO.ID,": ",CC_GO_Sub$Term)
MF_GO_Sub$GO_Term <- paste0(MF_GO_Sub$GO.ID,": ",MF_GO_Sub$Term)

#Combine results
GO_Sub <- rbind(CC_GO_Sub,MF_GO_Sub)

#Make Sankey diagram####

#Keep necessary columns for links
links <- GO_Sub[,c("Associated_DE_Gene","GO_Term")]

links$value <- 1

colnames(links)[1:2] <- c("target","source")

nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

nodes$group <- c("0","1","2","3","4","5","6","7","8","9",paste(rep("gene",48),sep=","))

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

links$group <- factor(links$IDsource)

#Set colors for links and nodes (cyan, gold, lavender, reddish pink, light orange,olive green, cornflower blue, bright pink, bright green, mint green,gray)
my_color <- 'd3.scaleOrdinal() .domain(["0","1","2","3","4","5","6","7","8","9","gene"]) .range(["#66cccc","#ffcc33","#996699","#ff6666","#f4896c","#cccc66","#99ccff","#ff66cc","#99ff00","#99ffcc","#999999"])'

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource",Target = "IDtarget", Value = "value", NodeID = "name", colourScale=my_color, NodeGroup="group", LinkGroup="group", nodeWidth = 20,fontSize = 16,  fontFamily = "sans-serif",sinksRight=FALSE,margin = list(top=4, right=2, bottom=4, left=2))

saveNetwork(p, "sankey_diagram_GO_Analysis2.html", selfcontained = TRUE)
