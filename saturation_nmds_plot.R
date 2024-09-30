# Code by: Isabelle Du Plessis, 2023 and Marian Dominguez-Mirazo, 2023

require(ggplot2) # version 3.4.2
require(reshape2) # version 1.4.4
require(vegan) # version 2.5-6
require(grid) # version 4.0.2
require(gridExtra) # version 2.3

# R version 4.0.2 (2020-06-22)

#### Load Data ####

path = getwd() # Set path

# Load sample metadata
metadata = read.csv(paste0(path, "/intermediate_files/metadata.csv"))
# Change sample id so that we can compare
metadata$Sample_ID = paste("sample_",metadata$Sample_ID,sep="")
# Change order so that we can compare
metadata = metadata[order(metadata$Sample_ID),]

# Load vOTU data
votus_cov75thres = read.delim(paste0(path,"/intermediate_files/votus_cov75thres.txt"), header=FALSE, comment.char="#")
colnames(votus_cov75thres) = c("sample","votus","coverage","meandepth","rpkm")

##############################################
## Cumulative vOTUs / Saturation Plot
perms = 100
le_tmp = table(votus_cov75thres[,1:2])

main = matrix(0,nrow=perms,ncol=24)
for(perm in 1:perms){
  tots = rep(0,769)
  tmp = le_tmp[sample(1:24),]
  for(i in 1:24){
    tmp2= tmp[i,] - tots
    tmp2[tmp2<0]=0
    main[perm,i] = sum(tmp2) #change the 1
    tots = tots + tmp2
    tots[tots>1]=1
  }
  main[perm,] = cumsum(main[perm,])
}
lemean = colSums(main)/perms
df_mean = data.frame("sample"=1:24,"means"=lemean)
df = melt(main)

saturation = ggplot(df_mean,aes(sample, means)) +
  geom_point(data=df, aes(Var2, value, group='Var1'),color = "darkgrey") +
  geom_line(color="dodgerblue2", linewidth = 1) +
  xlab("Number of Samples") +
  scale_y_continuous("Cumulative Number of vOTUs") +
  theme_test() + ggtitle("") + 
  theme(text = element_text(size = 15), plot.title = element_text(face = "bold", size = 15, hjust = .5),
        plot.margin = unit(c(5.5, 20, 5.5, 5.5),"pt"))



##############################################


### NMDS
sample = votus_cov75thres$sample
votu = votus_cov75thres$votus
rpkm = votus_cov75thres$rpkm
newdf2 = data.frame(sample,votu,rpkm)

votu = unique(newdf2$votu)
samples = sort(unique(newdf2$sample)) 
df1 = matrix(nrow = length(samples), ncol = length(votu), dimnames = list(samples,votu))


for(r in 1:nrow(newdf2)){
  samp = newdf2[r, 1]
  tax = newdf2[r, 2]
  df1[samp,tax] = newdf2[r, 3]
} # 1, 2, 3 here relate the the column number in the raw data in which the sample name, species name and data are in

df1[is.na(df1)] = 0   #convert NA's to 0
#df2 = as.data.frame(df1)

nmds = metaMDS(df1, distance = "bray", autotransform = FALSE)
#nmds #below 0.2 is generally good 
#plot(nmds) #circles are samples, starts are votus

data.scores = as.data.frame(scores(nmds))
#data.scores = as.data.frame((scores(nmds)$sites)) #problem line

data.scores$Sample = samples
height_list = c()
compartment_list = c()
for (i in (data.scores$Sample)){
  if (i %in% (c("sample_1","sample_7","sample_13","sample_19"))){
    compartment_list = append(compartment_list,"Bulk")
    height_list = append(height_list,"Tall")
  }
  if (i %in% (c("sample_2","sample_8","sample_14","sample_20"))){
    compartment_list = append(compartment_list,"Rhizosphere")
    height_list = append(height_list,"Tall")
  }
  if (i %in%(c("sample_3","sample_9","sample_15","sample_21"))){
    compartment_list = append(compartment_list, "Endosphere")
    height_list = append(height_list, "Tall")
  }
  if (i %in%c("sample_4","sample_10","sample_16","sample_22")){
    compartment_list = append(compartment_list, "Bulk")
    height_list = append(height_list, "Short")
  }
  if (i %in%(c("sample_5","sample_11","sample_17","sample_23"))){
    compartment_list = append(compartment_list, "Rhizosphere")
    height_list = append(height_list, "Short")
  }
  if (i %in% c("sample_6","sample_12","sample_18","sample_24")){
    compartment_list = append(compartment_list, "Endosphere")
    height_list = append(height_list, "Short")
  }
}
data.scores$Compartment = compartment_list
data.scores$Height = height_list

nmds_plot = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 2.5, aes( shape = Height, colour = Compartment)) + theme_test() +
  scale_colour_manual(values = c("Bulk" = "#853512", "Endosphere" = "#558A78","Rhizosphere" = "#EEAA23"), limits = c("Bulk", "Rhizosphere", "Endosphere"), labels = c("Bulk Sediment", "Rhizosphere", "Root")) +
  ggtitle("") + labs(color = "", shape = "") +
  theme(text = element_text(size = 15), plot.title = element_text(face = "bold", size = 15, hjust = .5))


# Plot

sat = grid.arrange(saturation, top = grid::textGrob("b", gp=gpar(fontsize=15,font=2), x = 0, hjust = 0))
nm = grid.arrange(nmds_plot, top = grid::textGrob("c", gp=gpar(fontsize=15,font=2), x = 0, hjust = 0))

grid.arrange(sat,nm, nrow=1, ncol=2, widths = c(2.5, 3.4))

