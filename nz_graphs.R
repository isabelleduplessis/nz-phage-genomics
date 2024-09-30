# New Zealand Viral Sequence Graphs
# Code by: Isabelle Du Plessis, 2023

require(ggplot2)
require(reshape2)
require(gplots)
require(vegan)
library(readr)
library(VennDiagram)
require(tidyr)
require(ggdendro)
require(ggvenn)
require(ggVennDiagram)



# sequence length distribution
nzseqlens <- read.table("~/phage/NZ/abundance/nzseqlens.txt", quote="\"", comment.char="")

# include all sequences
hist(nzseqlens$V1, breaks=max(nzseqlens$V1)/100)

# include only sequences smaller than 50,000 bp
hist(nzseqlens$V1[nzseqlens$V1<50000], breaks=50000/100, xlab = "Sequence Length (bin size = 100)", main = "Sequence Length Distribution") 
# still need to try to add cut ahead axis



#load sample metadata
virome_with_bottle <- read.csv("~/phage/NZ/abundance/virome_with_bottle.csv")
metadata <- data.frame("Sample"=virome_with_bottle$sample_id, "Depth"=virome_with_bottle$depth_numeric, "Cycle"=virome_with_bottle$Cycle, "Temperature"=virome_with_bottle$Temp1...C., "Salinity"=virome_with_bottle$Salinity1.PSS.78., "Oxygen"=virome_with_bottle$Oxygen1..mmol.kg., "Water Mass"=virome_with_bottle$water_mode)

metadata$Sample=parse_number(metadata$Sample)
metadata=metadata[metadata$Sample>=3,]
newmetadata=data.frame(aggregate(list(metadata$Depth, metadata$Cycle, metadata$Temperature, metadata$Salinity, metadata$Oxygen), by=list(metadata$Sample), FUN=mean))
newmetadata$Water="Subtropical"
colnames(newmetadata)=c("Sample", "Depth", "Cycle", "Avg Temperature", "Avg Salinity", "Avg Oxygen", "Water")
newmetadata$Water[newmetadata$Cycle!=3]="Subantarctic"


# proportion of reads that mapped to each sample
totalreads_mapped <- read.delim("~/phage/NZ/abundance/totalreads_mapped.csv", header=FALSE)
totalreads <- read.csv("~/phage/NZ/abundance/totalreads.csv", header=FALSE)
totalreads = totalreads[order(totalreads$V1),]
totalreads=totalreads[parse_number(totalreads$V1)>=3,]
viralreads = totalreads_mapped
viralreads = viralreads[order(viralreads$V1),]
viralreads=viralreads[parse_number(viralreads$V1)>=3,]

propviralreads = data.frame("Sample"=parse_number(viralreads$V1),"Viral"=viralreads$V2/totalreads$V2, "Water"=newmetadata$Water)
propviralreads$Water = factor(propviralreads$Water,levels = c("Subtropical", "Subantarctic"))
colnames(propviralreads) = c("Sample","Viral", "Water")

ggplot(data=propviralreads, aes(y =Viral, x=Water, color=Water)) +
  geom_point() +
  theme_test() +
  xlab("Sample") +
  ylab("Proportion of sample reads that map to viral sequence")



############################
## Total RPKM 

votus_cov75thres <- read.delim("~/phage/NZ/abundance/nz_votus_cov75thres_afterincreasing.txt", header=FALSE, comment.char="#")
#votus_cov75thres <- read.delim("~/Desktop/nz_votus_cov75thres_beforeincreasing.txt", header=FALSE, comment.char="#")
colnames(votus_cov75thres) = c("sample","votus","coverage","meandepth","rpkm")

#Plot total rpkm per sample, grouped by metadata
tmp = as.numeric(as.character(votus_cov75thres$rpkm))
tmp = aggregate(tmp,list(sample=votus_cov75thres$sample),sum)
colnames(tmp)=c("sample","rpkm")
tmp$sample=parse_number(tmp$sample)
tmp=tmp[tmp$sample>=3,]
tmp=tmp[order(tmp$sample),]
tmp = cbind(tmp, newmetadata[,2:7])
tmp$Water = factor(tmp$Water,levels = c("Subtropical", "Subantarctic"))
tmp$Cycle=as.character(tmp$Cycle)


rpkmpersample = tmp

ggplot(data=tmp,aes(Water,rpkm,color=Water)) +
  #geom_boxplot(color="black") +
  geom_point(size=2.5) +
  theme_test() + 
  xlab("Water Type") + 
  ylab("Total RPKM") +
  theme(axis.text.x=element_text(size=13), axis.text.y=element_text(size=13), legend.text = element_text(size=13), axis.title.x = element_text(size=15),axis.title.y = element_text(size=15))


# rpkm vs salinity colored by cycle
ggplot(data=tmp,aes(`Avg Salinity`,rpkm,color=Cycle)) +
  #geom_boxplot(color="black") +
  geom_point(size=2.5) +
  theme_test() + 
  xlab("Avg Salinity") + 
  ylab("Total RPKM") +
  theme(axis.text.x=element_text(size=13), axis.text.y=element_text(size=13), legend.text = element_text(size=13), axis.title.x = element_text(size=15),axis.title.y = element_text(size=15))

# rpkm vs temperature colored by cycle
ggplot(data=tmp,aes(`Avg Temperature`,rpkm,color=Cycle)) +
  #geom_boxplot(color="black") +
  geom_point(size=2.5) +
  theme_test() + 
  xlab("Avg Temperature") + 
  ylab("Total RPKM") +
  theme(axis.text.x=element_text(size=13), axis.text.y=element_text(size=13), legend.text = element_text(size=13), axis.title.x = element_text(size=15),axis.title.y = element_text(size=15))


############################
#Number of votus per sample
tmp = melt(table(votus_cov75thres$sample))
colnames(tmp) = c("sample","votus")
tmp$sample=parse_number(as.character(tmp$sample))
tmp=tmp[tmp$sample>=3,]
tmp=tmp[order(tmp$sample),]
tmp = cbind(tmp, newmetadata[,2:7])
tmp$Water = factor(tmp$Water,levels = c("Subtropical", "Subantarctic"))
tmp$Cycle=as.character(tmp$Cycle)



newtable = data.frame(rpkmpersample$sample, rpkmpersample$rpkm, tmp$votus, rpkmpersample$Depth, rpkmpersample$Cycle)
colnames(newtable) = c("Sample", "RPKM", "vOTUs", "Depth", "Cycle")




ggplot(data=tmp,aes(Water,votus,color=Water)) +
  geom_point(size=2.5) +
  theme_test() + 
  xlab("Water Type") + 
  ylab("Number of VOTUs Present")+
  theme(axis.text.x=element_text(size=13), axis.text.y=element_text(size=13), legend.text = element_text(size=13), axis.title.x = element_text(size=15),axis.title.y = element_text(size=15))

# votus vs salinity
ggplot(data=tmp,aes(`Avg Salinity`,votus,color=Cycle)) +
  geom_point(size=2.5) +
  theme_test() + 
  xlab("Avg Salinity") + 
  ylab("Number of VOTUs Present")+
  theme(axis.text.x=element_text(size=13), axis.text.y=element_text(size=13), legend.text = element_text(size=13), axis.title.x = element_text(size=15),axis.title.y = element_text(size=15))

# votu vs temperature
ggplot(data=tmp,aes(`Avg Temperature`,votus,color=Cycle)) +
  geom_point(size=2.5) +
  theme_test() + 
  xlab("Avg Temperature") + 
  ylab("Number of VOTUs Present")+
  theme(axis.text.x=element_text(size=13), axis.text.y=element_text(size=13), legend.text = element_text(size=13), axis.title.x = element_text(size=15),axis.title.y = element_text(size=15))



############################
#### Shannon diversity for full mapping
df = data.frame(matrix(ncol = 2, nrow = 20))
colnames(df) = c("sample","diversity")

votus_cov75thres=votus_cov75thres[parse_number(votus_cov75thres$sample)>=3,]

for(i in 3:22){
  sampleid = i
  for_shannon = votus_cov75thres[parse_number(votus_cov75thres$sample)==sampleid,]
  this_div  = diversity(for_shannon$rpkm)
  df[i,] = c(sampleid,this_div)
}

df=df[! is.na(df$sample),]
df = df[order(df$sample),]
df = cbind(df, newmetadata[,7])
df$diversity = as.numeric(as.character(df$diversity))
colnames(df)=c("sample", "diversity", "Water")
df$Water = factor(df$Water,levels = c("Subtropical", "Subantarctic"))

ggplot(data=df,aes(Water,diversity,color=Water)) +
  #geom_boxplot(color="black") +
  geom_point(size=2.5) +
  theme_test() + 
  xlab("") + 
  ylab("Shannon Diversity index")



## Cumulative vOTUs

perms = 100
le_tmp = table(votus_cov75thres[,1:2])
main = matrix(0,nrow=perms,ncol=20)
for(perm in 1:perms){
  tots = rep(0,7890)
  tmp = le_tmp[sample(1:20),]
  for(i in 1:20){
    tmp2= tmp[i,] - tots
    tmp2[tmp2<0]=0
    main[perm,i] = sum(tmp2) #change the 1
    tots = tots + tmp2
    tots[tots>1]=1
  }
  main[perm,] = cumsum(main[perm,])
}
lemean = colSums(main)/perms
df_mean = data.frame("sample"=1:20,"means"=lemean)
df = melt(main)
ggplot(df,aes(Var2,value,group='Var1')) +
  geom_point() +
  geom_line(data=df_mean,aes(sample,means),color="blue") +
  xlab("Number of Samples") +
  ylab("Cumulative Number of vOTUs") +
  theme_test()



# Rank order graph
my_table = votus_cov75thres[,c(1,2,5)]
for(i in 3:22){
  my_sample = i
  my_table2 = my_table[parse_number(my_table$sample) == my_sample,]
  my_table2 = my_table2[order(my_table2$rpkm,decreasing = T),]
  nam <- paste("df", i, sep = "")
  tmpdf = data.frame(1:nrow(my_table2), cumsum(my_table2$rpkm))
  colnames(tmpdf) = c("order","rpkm")
  tmpdf$sample=rep(i, nrow(tmpdf))
  tmpdf$rpkm=tmpdf$rpkm/max(tmpdf$rpkm)
  assign(nam, tmpdf)
}


totaldf = rbind(df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14, df15, df16, df17, df18, df19, df20, df21, df22)
totaldf$sample = as.character(totaldf$sample)
totaldf$order=totaldf$order/max(totaldf$order) # normalize x axis zero to 1 percentage

ggplot(totaldf,aes(order,rpkm,color=sample)) + 
  geom_line() +
  theme_test() 



#############################
# NMDS plot
votus_cov75thres=votus_cov75thres[parse_number(votus_cov75thres$sample)>=3,]
otherdf=votus_cov75thres[,c(1,2,5)]
rpkm = dcast(otherdf, sample ~ votus)
rpkm[is.na(rpkm)] <- 0
rpkm$sample=parse_number(rpkm$sample)

rpkm=rpkm[rpkm$sample!=19]

m_rpkm = as.matrix(rpkm)

k5 = kmeans(m_rpkm, 5, nstart=25)
fviz_cluster(k5, m_rpkm, ggtheme = theme_test(), geom = "point")

nmds = metaMDS(m_rpkm, distance="bray")
data.scores = as.data.frame(scores(nmds))
data.scores$Sample = rpkm$sample
data.scores=data.scores[order(data.scores$Sample),]
tmp = as.numeric(as.character(votus_cov75thres$rpkm[votus_cov75thres$sample!=19]))
#tmp = as.numeric(as.character(votus_cov75thres$rpkm))
tmp = aggregate(tmp,list(sample=votus_cov75thres$sample[votus_cov75thres$sample!=19]),sum)
#tmp = aggregate(tmp,list(sample=votus_cov75thres$sample),sum)
colnames(tmp)=c("sample","rpkm")
tmp$sample=parse_number(tmp$sample)
tmp=tmp[tmp$sample>=3,]
tmp=tmp[order(tmp$sample),]
#newmetadata=newmetadata[newmetadata$sample!=19,]
tmp = cbind(tmp, newmetadata[ ,2:7])
tmp$Water = factor(tmp$Water,levels = c("Subtropical", "Subantarctic"))
data.scores=data.scores[data.scores$Sample>=3,]
data.scores$Water=tmp$Water
data.scores$Cycle=tmp$Cycle

data.scores=data.scores[data.scores$Sample!=19,]

ggplot(data.scores, aes(x=NMDS1, y=NMDS2, shape = as.character(Cycle), color = Water)) +
  geom_point(size = 4) + 
  labs(x = "NMDS1", shape = "Cycle", y = "NMDS2", colour = "Water") +
  theme_test()



# cluster by presence or absence
votus_cov75thres=votus_cov75thres[parse_number(votus_cov75thres$sample)>=3,]
otherdf=votus_cov75thres[,c(1,2)]
otherdf$value = 1
rpkm = dcast(otherdf, sample ~ votus)
rpkm[is.na(rpkm)] <- 0
rpkm$sample=parse_number(rpkm$sample)





library(ggpubr)
library(factoextra)

m_rpkm = as.matrix(rpkm)

k5 = kmeans(m_rpkm, 5, nstart=25)
fviz_cluster(k5, m_rpkm, ggtheme = theme_test(), geom = "point")

## heatmap - includes sample 19
votus_cov75thres=votus_cov75thres[parse_number(votus_cov75thres$sample)>=3,]
otherdf=votus_cov75thres[,c(1,2,5)]
rpkm = dcast(otherdf, sample ~ votus) 
rpkm[is.na(rpkm)] <- 0
rpkm$sample=parse_number(rpkm$sample)
rpkm = rpkm[order(rpkm$sample),]

rpkm=rpkm[,1:150]
# for each sample, sum all the viruses in the sample and divide each by total number to normalize
for(i in 1:20){
  #totalsum = sum(rpkm[i,-1])
  rpkm[i,-1]=rpkm[i,-1]/max(rpkm[i,-1])
}



h = hclust(dist(rpkm[,-1]))
h$order
plot(h)

ord = c(h$order)

dendro <- as.dendrogram(h)

# Create dendro
dendro_plot <- ggdendrogram(data = dendro, rotate = TRUE)

# Preview the plot
print(dendro_plot)

long <- pivot_longer(data = rpkm, cols=-c(sample))

long = long[long$sample!=19,]


ggplot(data = long, aes(x = name, y = as.character(sample))) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2() +
  theme(axis.text.y = element_text(size = 6)) +
  scale_y_discrete(limits=c("3", "4", "19", "20", "17", "8", "18", "5", "6", "7", "11", "12", "16", "9", "10", "15", "13", "14")) +
  scale_x_discrete(labels=NULL) +
  xlab("vOTUs") +
  ylab("Sample")


## venn diagram

newdf = votus_cov75thres[,c(1,2)]


c1 = votus_cov75thres[parse_number(votus_cov75thres$sample) %in% c(3, 4, 5, 6),] #samples 3, 4, 5, 6
c2 = votus_cov75thres[parse_number(votus_cov75thres$sample) %in% c(7, 8, 9, 10),] #samples 7, 8, 9, 10
c3 = votus_cov75thres[parse_number(votus_cov75thres$sample) %in% c(11, 12, 13, 14),] #samples 11, 12, 13, 14
c4 = votus_cov75thres[parse_number(votus_cov75thres$sample) %in% c(15, 16, 17, 18),] #samples 15, 16, 17, 18
c5 = votus_cov75thres[parse_number(votus_cov75thres$sample) %in% c(19, 20, 21, 22),] #samples 19, 20, 21, 22



###
c1$cycle = as.character(1)
c2$cycle = as.character(2)
c3$cycle = as.character(3)
c4$cycle = as.character(4)
c5$cycle = as.character(5)

totalcycles = unique(rbind(c1[,c(2,6)], c2[,c(2,6)], c3[,c(2,6)], c4[,c(2,6)], c5[,c(2,6)]))
totalcycles$val = 1


plotcycles = dcast(totalcycles, votus ~ cycle)

plotcycles[is.na(plotcycles)] <- 0

colnames(plotcycles) = c("votus", "cycle1", "cycle2", "cycle3", "cycle4", "cycle5")


ggplot(plotcycles, aes(A = as.logical(cycle1), B = as.logical(cycle2), C = as.logical(cycle3), D = as.logical(cycle4), E = as.logical((cycle5)))) +
  geom_venn(show_percentage = FALSE, set_names = c("Cycle 1", "Cycle 2", "Cycle 3", "Cycle 4", "Cycle 5")) + 
  theme_test() +
  theme(text = element_text(size = 15), panel.border = element_blank()) +
  scale_x_continuous(NULL, NULL, NULL) +
  scale_y_continuous(NULL, NULL, NULL)

####








library(RColorBrewer)
myCol <- brewer.pal(9, "YlGnBu")




myCol = c("#41B6C4", "#FC8D62", "#A6D854", "#E78AC3", "#FFD92F")


venn.diagram(x=list(unique(c1$votus), unique(c2$votus), unique(c3$votus), unique(c4$votus), unique(c5$votus)), 
             category.names = c("Cycle 1", "Cycle 2", "Cycle 3", "Cycle 4", "Cycle 5"), 
             filename = 'cycles.png',
             output=TRUE,
             #lty = 'blank',
             fill = myCol,
             width = 4000, 
             cat.fontface = "bold",
             lwd = 2,
             fontfamily = "sans",
             cat.fontfamily = "sans",
             cat.pos = c(0,0,0,0,0),
             cat.dist = c(.2, .2, -.15, -.2,.2),
             cex = .9)

