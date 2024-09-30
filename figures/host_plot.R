# Code by: Isabelle Du Plessis, 2023

require(ggplot2) # version 3.4.2
require(readr) # version 2.1.4
require(gridExtra) # version 2.3
require(ggpubr) # version 0.6.0
require(ggvenn) # version 0.1.10

# R version 4.0.2 (2020-06-22)

#### Load Data ####

path = getwd() # Set path

# Load sample metadata
metadata = read.csv(paste0(path, "/intermediate_files/metadata.csv"))

# Load vOTU data
votus_cov75thres = read.delim(paste0(path,"/intermediate_files/votus_cov75thres.txt"), header=FALSE, comment.char="#")
colnames(votus_cov75thres) = c("sample","votus","coverage","meandepth","rpkm")
votus_cov75thres$sample = parse_number(votus_cov75thres$sample)

# Load host data
HostPrediction = read.delim(paste0(path, "/intermediate_files/HostPrediction.tsv"))
hostdf = read.table(paste0(path, "/intermediate_files/hostproperties.tsv"), quote="\"", comment.char="")
colnames(hostdf) = c("virus", "hostorder", "sulfuroxidizing", "sulfatereducing", "ironoxidizing", "nitrifying")


#### Create main dataframe ####

# Create dataframe with all vOTUs and all info
bigdf = data.frame(votus_cov75thres[, c(1,2)])
colnames(bigdf) = c("sample","virus")

# Initialize host vectors with domain, phylum, and order values
d <- rep("Unknown", nrow(bigdf))
p <- rep("Unknown", nrow(bigdf))
o <- rep("Unknown", nrow(bigdf))

# For each vOTU, find corresponding compartment, phenotype, and host prediction
for(i in 1:nrow(bigdf)){
  bigdf$Compartment[i]=metadata$Compartment[bigdf$sample[i]==metadata$Sample_ID] # Add compartment to bigdf
  bigdf$Phenotype[i]=metadata$Spartina[bigdf$sample[i]==metadata$Sample_ID] # Add phenotype to bigdf
  if(bigdf$virus[i] %in% HostPrediction$Virus){ # For vOTUs that have hosts predicted, add host taxonomy to host vectors
    index=which(bigdf$virus[i]==HostPrediction$Virus)
    d[i]=HostPrediction$Domain[index]
    p[i]=HostPrediction$Phylum[index]
    o[i]=HostPrediction$Order[index]
  }
}

# Set new bigdf columns equal to host vectors
bigdf$Domain=d 
bigdf$Phylum=p
bigdf$Order=o

# Factor compartments for plotting
bigdf$Compartment = factor(bigdf$Compartment,levels = c("Bulk sediment","Rhizosphere","Endosphere"))

# Initialize property columns
bigdf$sulfuroxidizing = 0
bigdf$sulfatereducing = 0
bigdf$ironoxidizing = 0

# Find corresponding biogeochemical properties
for(i in 1:nrow(bigdf)){ 
  virus = bigdf[i,2] 
  if(virus %in% hostdf$virus){ # For each virus in bigdf that has a predicted host, add host biogeochemical properties
    bigdf[i,8] = hostdf[which(hostdf$virus == virus),3]
    bigdf[i,9] = hostdf[which(hostdf$virus == virus),4]
    bigdf[i,10] = hostdf[which(hostdf$virus == virus),5]
  }
}

# Main dataframe is complete

#### Phylum Stacked Bar Plot - Figure 5a ####

# Get phyla of predicted hosts, removing replicate viruses in the same compartment/phenotype
phylumdf = bigdf[,c(2,3,4,6)] # Get relevant columns: Compartment, Phenotype, and Phylum
phylumdf = unique(phylumdf[which(phylumdf$Phylum!="Unknown"),]) # Get unique rows to remove replicates, get only known hosts
phylumdf$value = 1 # For plotting

# Get proportion of vOTUs with hosts predicted out of all vOTUs per compartment/phenotype
phylumdf$hostcount = 0 # Initialize column for number of vOTUs with hosts predicted

# Get number of vOTUs with each phylum host prediction, per compartment/phenotype
phylumdf[which(phylumdf$Compartment=="Bulk sediment" & phylumdf$Phenotype=="Tall"),6] = nrow(phylumdf[which(phylumdf$Compartment=="Bulk sediment" & phylumdf$Phenotype=="Tall"),])
phylumdf[which(phylumdf$Compartment=="Bulk sediment" & phylumdf$Phenotype=="Short"),6] = nrow(phylumdf[which(phylumdf$Compartment=="Bulk sediment" & phylumdf$Phenotype=="Short"),])
phylumdf[which(phylumdf$Compartment=="Rhizosphere" & phylumdf$Phenotype=="Tall"),6] = nrow(phylumdf[which(phylumdf$Compartment=="Rhizosphere" & phylumdf$Phenotype=="Tall"),])
phylumdf[which(phylumdf$Compartment=="Rhizosphere" & phylumdf$Phenotype=="Short"),6] = nrow(phylumdf[which(phylumdf$Compartment=="Rhizosphere" & phylumdf$Phenotype=="Short"),])
phylumdf[which(phylumdf$Compartment=="Endosphere" & phylumdf$Phenotype=="Tall"),6] = nrow(phylumdf[which(phylumdf$Compartment=="Endosphere" & phylumdf$Phenotype=="Tall"),])
phylumdf[which(phylumdf$Compartment=="Endosphere" & phylumdf$Phenotype=="Short"),6] = nrow(phylumdf[which(phylumdf$Compartment=="Endosphere" & phylumdf$Phenotype=="Short"),])

phylumdf$totalcount = 0 # Initialize column for number of vOTUs total
tmp = unique(bigdf[,c(2,3,4)]) # Get unique viruses, compartment, and phenotype

# Get total number of vOTUs in each compartment/phenotype
phylumdf[which(phylumdf$Compartment=="Bulk sediment" & phylumdf$Phenotype=="Tall"),7] = nrow(tmp[which(tmp$Compartment=="Bulk sediment" & tmp$Phenotype=="Tall"),])
phylumdf[which(phylumdf$Compartment=="Bulk sediment" & phylumdf$Phenotype=="Short"),7] = nrow(tmp[which(tmp$Compartment=="Bulk sediment" & tmp$Phenotype=="Short"),])
phylumdf[which(phylumdf$Compartment=="Rhizosphere" & phylumdf$Phenotype=="Tall"),7] = nrow(tmp[which(tmp$Compartment=="Rhizosphere" & tmp$Phenotype=="Tall"),])
phylumdf[which(phylumdf$Compartment=="Rhizosphere" & phylumdf$Phenotype=="Short"),7] = nrow(tmp[which(tmp$Compartment=="Rhizosphere" & tmp$Phenotype=="Short"),])
phylumdf[which(phylumdf$Compartment=="Endosphere" & phylumdf$Phenotype=="Tall"),7] = nrow(tmp[which(tmp$Compartment=="Endosphere" & tmp$Phenotype=="Tall"),])
phylumdf[which(phylumdf$Compartment=="Endosphere" & phylumdf$Phenotype=="Short"),7] = nrow(tmp[which(tmp$Compartment=="Endosphere" & tmp$Phenotype=="Short"),])

# Plot phylum vs phenotype/compartment
phyla = ggplot(phylumdf, aes(interaction(Phenotype, Compartment), value, fill=Phylum)) +
  geom_bar(position = "fill", stat = "identity") +
  theme_test() +
  ylab("Proportion of Phylum-Level Hosts") +
  scale_fill_manual(values = c("#89C5DA", "#DA5724", "#CE50CA", "#8569D5", "#C0717C", 
                               "#689030", "#673770", "#5F7FC7", "#38333E", "#508578", "#C84248", "#CD9BCD", "#7FDCC0",
                               "#D1A33D")) +
  scale_x_discrete("Sample Compartment and Phenotype", labels=c("Bulk Sediment\nShort", "Bulk Sediment\nTall", "Rhizosphere\nShort", "Rhizosphere\nTall", "Root\nShort", "Root\nTall")) +
  theme(text = element_text(size = 15)) + theme(axis.text.x = element_text(size = 10)) +
  ggtitle("Predicted Target Host Phyla") + theme(plot.title = element_text(face = "bold", size = 15, hjust = .5)) +
  geom_text(aes(label = paste(round((hostcount/totalcount)*100), "%")), vjust = -0.4)


#### Host Property Venn Diagram - Figure 5b ####

# Plot venn diagram of properties shared by predicted hosts
host_venn = ggplot(hostdf, aes(A=as.logical(sulfuroxidizing), B=as.logical(sulfatereducing), C=as.logical(ironoxidizing))) + 
  geom_venn(show_percentage = FALSE, set_names = c("Sulfur Oxidizing  ", "    Sulfate Reducing", "Iron Oxidizing"), 
            fill_color = c("turquoise3", "brown", "gold"), text_size = 6, fill_alpha = .6) +
  theme_test() +
  theme(text = element_text(size = 15), panel.border = element_blank()) +
  scale_x_continuous(NULL, NULL, NULL) +
  scale_y_continuous(NULL, NULL, NULL) +
  ggtitle("Biogeochemical Properties of Host Orders") + theme(plot.title = element_text(face = "bold", size = 15, hjust = .5))

# Plot venn diagram of properties shared by unique predicted hosts
unique_host_venn = ggplot(unique(hostdf[,c(2,3,4,5)]), aes(A=as.logical(sulfuroxidizing), B=as.logical(sulfatereducing), C=as.logical(ironoxidizing))) + 
  geom_venn(show_percentage = FALSE, set_names = c("Sulfur Oxidizing  ", "    Sulfate Reducing", "Iron Oxidizing"), 
            fill_color = c("turquoise3", "brown", "gold"), text_size = 6, fill_alpha = .6) +
  theme_test() +
  theme(text = element_text(size = 15), panel.border = element_blank()) +
  scale_x_continuous(NULL, NULL, NULL) +
  scale_y_continuous(NULL, NULL, NULL) +
  ggtitle("Biogeochemical Properties of Unique Host Orders") + theme(plot.title = element_text(face = "bold", size = 15, hjust = .5))


#### Host Property Bar Plot - Figure 5c ####

# Create dataframe for proportions of hosts with properties for each compartment/phenotype
propdf = unique(bigdf[,c(3,4)])

# Initialize cols
propdf$numso = 0
propdf$numsr = 0
propdf$numio = 0
propdf$totalhosts = 0
propdf$soprop = 0
propdf$srprop = 0
propdf$ioprop = 0


for(i in 1:nrow(propdf)){ # For each compartment/phenotype combination, get proportion of hosts with properties out of predicted hosts
  tothosts = nrow(bigdf[which(bigdf[,3]==propdf[i,1] & bigdf[,4]==propdf[i,2]),]) # Total number of hosts predicted
  num_so = nrow(bigdf[which((bigdf[,3]==propdf[i,1] & bigdf[,4]==propdf[i,2]) & bigdf[,8]==1),]) # Number sulfur oxidizing
  num_sr = nrow(bigdf[which(bigdf[,3]==propdf[i,1] & bigdf[,4]==propdf[i,2] & bigdf[,9]==1),]) # Number sulfate reducing
  num_io = nrow(bigdf[which(bigdf[,3]==propdf[i,1] & bigdf[,4]==propdf[i,2] & bigdf[,10]==1),]) # Number iron oxidizing
  propdf[i,3] = num_so  
  propdf[i,4] = num_sr 
  propdf[i,5] = num_io  
  propdf[i,6] = tothosts
  propdf[i,7] = num_so/tothosts # Proportion sulfur oxidizing
  propdf[i,8] = num_sr/tothosts # Proportion sulfate reducing
  propdf[i,9] = num_io/tothosts # Proportion iron oxidizing
}

# Plot proportion of hosts with sulfur oxidizing properties
sohosts = ggplot(propdf, aes(interaction(Compartment, Phenotype), soprop, fill = Compartment)) +
  geom_bar(stat = "identity") +
  theme_test() +
  scale_fill_manual(values = c("Bulk sediment" = "#853512",
                               "Endosphere" = "#558A78",
                               "Rhizosphere" = "#EEAA23")) +
  scale_x_discrete(NULL, limits = c("Bulk sediment.Short", "Bulk sediment.Tall", "Rhizosphere.Short", "Rhizosphere.Tall", "Endosphere.Short", "Endosphere.Tall"),
                   labels = c("Short", "Tall","Short", "Tall","Short", "Tall")) +
  ylim(0,.25) + ylab("") +
  theme(text = element_text(size = 15)) +
  theme(plot.margin = unit(c(5.5, 0, 5.5, 0),"pt"))

# Plot proportion of hosts with sulfate reducing properties
srhosts = ggplot(propdf, aes(interaction(Compartment, Phenotype), srprop, fill = Compartment)) +
  geom_bar(stat = "identity") +
  theme_test() +
  scale_fill_manual(values = c("Bulk sediment" = "#853512",
                               "Endosphere" = "#558A78",
                               "Rhizosphere" = "#EEAA23")) +
  scale_x_discrete(NULL, limits = c("Bulk sediment.Short", "Bulk sediment.Tall", "Rhizosphere.Short", "Rhizosphere.Tall", "Endosphere.Short", "Endosphere.Tall"),
                   labels = c("Short", "Tall","Short", "Tall","Short", "Tall")) +
  ylim(0,.25) + ylab("") +
  theme(text = element_text(size = 15)) +
  theme(plot.margin = unit(c(5.5, 5.5, 5.5, 0),"pt"))

# Plot proportion of hosts with iron oxidizing properties
iohosts = ggplot(propdf, aes(interaction(Compartment, Phenotype), ioprop, fill = Compartment)) +
  geom_bar(stat = "identity") +
  theme_test() +
  scale_fill_manual(values = c("Bulk sediment" = "#853512",
                               "Endosphere" = "#558A78",
                               "Rhizosphere" = "#EEAA23"), limits = c("Bulk sediment", "Rhizosphere", "Endosphere"), labels = c("Bulk Sediment", "Rhizosphere", "Root")) +
  scale_x_discrete(NULL, limits = c("Bulk sediment.Short", "Bulk sediment.Tall", "Rhizosphere.Short", "Rhizosphere.Tall", "Endosphere.Short", "Endosphere.Tall"),
                   labels = c("Short", "Tall","Short", "Tall","Short", "Tall")) +
  ylim(0,.25) +
  ylab("Proportion of Order-Level Hosts") +
  theme(text = element_text(size = 15)) +
  theme(legend.position = "bottom", legend.box = "horizontal") + labs(fill = NULL) +
  theme(plot.margin = unit(c(5.5, 0, 5.5, 5.5),"pt"))



## Plot all

# Add labels to figure 5 a and b
fig5a = grid.arrange(phyla, top = textGrob("a", gp=gpar(fontsize=15,font=2), x = 0, hjust = 0))
fig5b = grid.arrange(unique_host_venn, top = textGrob("b", gp=gpar(fontsize=15,font=2), x = 0, hjust = 0))

# Arrange figure 5 c
p1 = grid.arrange((iohosts + theme(legend.position = "none")), top = textGrob("Iron Oxidizing", gp=gpar(fontsize=15,font=2)))
p2 = grid.arrange((sohosts + theme(legend.position = "none")), top = textGrob("Sulfur Oxidizing", gp=gpar(fontsize=15,font=2)))
p3 = grid.arrange((srhosts + theme(legend.position = "none")), top = textGrob("Sulfate Reducing", gp=gpar(fontsize=15,font=2)))

# Add labels to figure 5 c
fig5c = grid.arrange(p1, p2, p3, get_legend(iohosts), nrow=2,
                     layout_matrix= rbind(c(1,2,3), c(4,4,4)),
                     widths = c(2, 2, 2), heights = c(4,1), top = textGrob("Biogeochemical Properties of Host Orders by Compartment and Phenotype", gp=gpar(fontsize=15,font=2)))

fig5c = grid.arrange(fig5c, top = textGrob("c", gp=gpar(fontsize=15,font=2), x = 0, hjust = 0))

# Put figure 5 together
grid.arrange(fig5a, fig5b, fig5c, nrow = 2, ncol = 2, widths = c(3.7, 2), heights = c(1.2, 1), layout_matrix = rbind(c(1, 2), c(3, 3)))

