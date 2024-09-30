# Code by: Isabelle Du Plessis, 2023

require(ggplot2) # version 3.4.2
require(reshape2) # version 1.4.4
require(readr) # version 2.1.4
require(grid) # version 4.0.2
require(gridExtra) # version 2.3
require(ggpubr) # version 0.6.0

# R version 4.0.2 (2020-06-22)

#### Load Data
###############

path = getwd() # Set path

# load sample metadata
metadata = read.csv(paste0(path, "/intermediate_files/metadata.csv"))
#change sample id so that we can compare
metadata$Sample_ID = paste("sample_",metadata$Sample_ID,sep="")
#change order so that we can compare
metadata = metadata[order(metadata$Sample_ID),]

# load votu data
votus_cov75thres = read.delim(paste0(path,"/intermediate_files/votus_cov75thres.txt"), header=FALSE, comment.char="#")
colnames(votus_cov75thres) = c("sample","votus","coverage","meandepth","rpkm")

# load qPCR data
qPCR_data = read.csv(paste0(path, "/intermediate_files/qPCR_data.csv"))


##############################################
## Viral Abundance, total RPKM per sample, grouped by metadata

tmp = as.numeric(as.character(votus_cov75thres$rpkm))
tmp = aggregate(tmp,list(sample=votus_cov75thres$sample),sum)
colnames(tmp)=c("sample","rpkm")
tmp = cbind(tmp, metadata[,2:5])
tmp$Compartment = factor(tmp$Compartment,levels = c("Bulk sediment","Rhizosphere","Endosphere"))

viral_abundance = ggplot(data=tmp,aes(Compartment,rpkm,color=Compartment,group=Spartina)) +
  geom_point(aes(shape = Spartina), size=2.5, position = position_dodge(width = .9)) +
  theme_test() + 
  stat_summary(fun = mean, geom = "tile", color = "black", height = ((max(tmp$rpkm)-min(tmp$rpkm))*.003), width = .8, position = position_dodge(width = .9)) +
  scale_color_manual(values = c("Bulk sediment" = "#853512",
                                "Endosphere" = "#558A78",
                                "Rhizosphere" = "#EEAA23"), labels = c("Bulk Sediment", "Rhizsophere", "Root")) +
  scale_x_discrete(NULL, labels = c("Bulk Sediment", "Rhizsophere", "Root")) +
  ylab("Total Viral RPKM") +
  theme(text = element_text(size = 15)) +
  theme(legend.position = "none") +
  theme(plot.margin = unit(c(15, 15, 30, 15),"pt")) +
  ggtitle("Viral Abundance") + theme(plot.title = element_text(face = "bold", size = 15, hjust = .5))


##############################################
#Number of votus per sample

tmp = melt(table(votus_cov75thres$sample))
colnames(tmp) = c("sample","votus")
tmp = cbind(tmp, metadata[,2:5])
tmp$Compartment = factor(tmp$Compartment,levels = c("Bulk sediment","Rhizosphere","Endosphere"))


votu_graph = ggplot(data=tmp,aes(Compartment, votus, color=Compartment, group=Spartina)) +
  geom_point(aes(shape = Spartina), size=2.5, position = position_dodge(width = .9)) +
  theme_test() + 
  stat_summary(fun = mean, geom = "tile", color = "black", height = ((max(tmp$votus)-min(tmp$votus))*.003), width = .8, position = position_dodge(width = .9)) +
  scale_color_manual(values = c("Bulk sediment" = "#853512",
                                "Endosphere" = "#558A78",
                                "Rhizosphere" = "#EEAA23"), labels = c("Bulk Sediment", "Rhizosphere", "Root")) +
  scale_x_discrete(NULL, labels = c("Bulk Sediment", "Rhizosphere", "Root")) + 
  ylab("Number of vOTUs") +
  theme(text = element_text(size = 15)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 5.5, b = 0, l = 0))) +
  theme(legend.position = "left") +  theme(legend.box = "horizontal") +
  labs(color = NULL, shape = NULL) +
  theme(plot.margin = unit(c(15, 15, 30, 5.5),"pt")) +
  ggtitle("vOTU Count") + theme(plot.title = element_text(face = "bold", size = 15, hjust = .5))


##### Prokaryotic Abundance ##################
abundance = data.frame("Compartment" = qPCR_data$Compartment, "Phenotype" = qPCR_data$Spartina, "Log Copies" = qPCR_data$logcopies_g_mass)
abundance$Compartment[abundance$Compartment=="Root"]="Endosphere"
abundance$Compartment = factor(abundance$Compartment, levels = c("Bulk sediment","Rhizosphere","Endosphere"))


prok_abundance = ggplot(abundance, aes(x=Compartment, y=Log.Copies, color = Compartment, group = Phenotype)) +
  geom_point(aes(shape = Phenotype), size=2.5, position = position_dodge(width = .9)) +
  theme_test() +
  scale_x_discrete(NULL, labels = c("Bulk Sediment", "Rhizosphere", "Root")) +
  scale_color_manual(values = c("Bulk sediment" = "#853512",
                                "Endosphere" = "#558A78",
                                "Rhizosphere" = "#EEAA23"), labels = c("Bulk Sediment", "Rhizsophere", "Root")) +
  ylab("Log(copies 16S rRNA gene) g-1") +
  stat_summary(fun = mean, geom = "tile", color = "black", height = ((max(abundance$Log.Copies)-min(abundance$Log.Copies))*.003), width = .8, position = position_dodge(width = .9)) +
  theme(text = element_text(size = 15)) + theme(legend.position = "none") +
  theme(plot.margin = unit(c(5.5, 15, 30, 15),"pt")) +
  ggtitle("Prokaryotic Abundance") + theme(plot.title = element_text(face = "bold", size = 15, hjust = .5))




##############################################
# Rank order 

my_table = votus_cov75thres[,c(1,2,5)]
my_table$Compartment = 0
my_table$Phenotype = 0

# for each votu, find compartment and phenotype based on sample number
for(i in 1:nrow(my_table)){
  my_table$Compartment[i]=metadata$Compartment[parse_number(my_table$sample[i])==parse_number(metadata$Sample_ID)]
  my_table$Phenotype[i]=metadata$Spartina[parse_number(my_table$sample[i])==parse_number(metadata$Sample_ID)]
  if(my_table$Compartment[i]=="Bulk sediment"){
    my_table$Compartment[i]="Bulk"
  }
}


my_table = my_table[,c(2,3,4,5)]


rowzero = data.frame(matrix(ncol = 4, nrow = 1))
colnames(rowzero) = c("order", "rpkm", "Compartment", "Phenotype")

# create a dataframe for each Compartment/Phenotype combination
for(i in c("Bulk", "Rhizosphere", "Endosphere")){
  for(j in c("Tall", "Short")){
    nam <- paste(i,j, sep="")
    
    
    
    my_table2 = my_table[which(my_table$Compartment == i),]
    my_table2 = my_table2[which(my_table2$Phenotype==j),]
    
    # Set rpkm to the sum of each set of replicates
    for(votu in unique(my_table2$votus)){
      indices = which(my_table2$votus == votu)
      tmp_sum = sum(my_table2[indices,2])
      my_table2[which(my_table2$votus == votu),]$rpkm = tmp_sum
    }

    # Get unique rows
    my_table2 = unique(my_table2)
    
    my_table2 = my_table2[order(my_table2$rpkm,decreasing = T),] # order rpkms
    tmpdf = data.frame(1:nrow(my_table2), cumsum(my_table2$rpkm)) # get cumulative sum of rpkms and assign order number
    colnames(tmpdf) = c("order","rpkm")
    tmpdf$Compartment=rep(i, nrow(tmpdf)) # add compartment column
    tmpdf$Phenotype=rep(j, nrow(tmpdf)) # add phenotype column
    tmpdf$rpkm=tmpdf$rpkm/max(tmpdf$rpkm) # make RPKM [0,1]
    for(k in 1:nrow(tmpdf)){
      tmpdf$order[k]=tmpdf$order[k]/max(tmpdf$order) # make order [0,1]
    }
    
    rowzero$order = 0
    rowzero$rpkm = 0
    rowzero$Compartment = i
    rowzero$Phenotype = j
    tmpdf2 = rbind(rowzero,tmpdf)
    assign(nam, tmpdf2)
  }
}

# combine df for each combination for plotting
totaldf = rbind(BulkTall, BulkShort, RhizosphereShort, RhizosphereTall, EndosphereShort, EndosphereTall)

rank_order_plot = ggplot(totaldf,aes(order,rpkm,group = interaction(Compartment, Phenotype), color = Compartment, linetype = Phenotype, label = )) + 
  geom_line(linewidth = .7) +
  scale_linetype_manual(values = c("Tall" = "solid", "Short" = "dashed")) +
  theme_test() +
  scale_color_manual(values = c("Bulk" = "#853512",
                                "Endosphere" = "#558A78",
                                "Rhizosphere" = "#EEAA23"), limits = c("Bulk", "Rhizosphere", "Endosphere"), labels = c("Bulk Sediment", "Rhizosphere", "Root")) +
  ylab("Relative Cumulative Abundance") + xlab("Proportion vOTU Rank") +
  theme(text = element_text(size = 15)) +
  theme(legend.box = "horizontal") + labs(color = NULL, linetype = NULL) +
  theme(legend.position = "left") + theme(plot.margin = unit(c(7.5, 5.5, 12, 5.5),"pt")) +
  ggtitle("vOTU Rank Order") + theme(plot.title = element_text(face = "bold", size = 15, hjust = .5))


########### plot #############


fig3a = grid.arrange(votu_graph + theme(legend.position = "none"), top = grid::textGrob("a", gp=gpar(fontsize=12,font=2), x = 0, hjust = 0))

fig3b = grid.arrange(viral_abundance + theme(legend.position = "none"), top = grid::textGrob("b", gp=gpar(fontsize=12,font=2), x = 0, hjust = 0))


fig3c = grid.arrange((rank_order_plot + theme(legend.position = "none")), top = textGrob("c", gp=gpar(fontsize=12,font=2), x = 0, hjust = 0))
fig3d = grid.arrange(prok_abundance, top = textGrob("d", gp=gpar(fontsize=12,font=2), x = 0, hjust = 0))

fig3cd <- grid.arrange(fig3c, fig3d, nrow = 1)

dot_legend <- get_legend(votu_graph)
line_legend <- get_legend(rank_order_plot)


grid.arrange(fig3a, fig3b, fig3cd, dot_legend, line_legend,
             layout_matrix= rbind(c(1,2), c(3, 3), c(5,4)), 
             heights = c(3,3,1), widths = c(1,1))

