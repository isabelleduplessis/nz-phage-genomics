# New Zealand Viral Sequence Graphs
# Code by: Isabelle Du Plessis, 2024 and Marian Dominguez-Mirazo, 2023
# R version 4.3.2

require(ggplot2) # v3.4.4
require(reshape2) # v1.4.4
require(readr) # v2.1.4
require(vegan) # v2.6.8
require(factoextra) # v1.0.7
require(dplyr) # v1.1.4
require(ggvenn) # v0.1.10


### Load Data ###

# Get path
path = getwd() # Set path

# Load sequence lengths
nzseqlens = read.table(paste0(path, "/data/nzseqlens.txt"), quote="\"", comment.char="")

# Load sample metadata
metadata = read.csv(paste0(path, "/data/metadata.csv"))

# Format metadata
metadata = data.frame("Sample"=metadata$sample_id, "Depth"=metadata$depth_numeric, "Cycle"=metadata$Cycle, "Temperature"=metadata$Temp1...C., "Salinity"=metadata$Salinity1.PSS.78., "Oxygen"=metadata$Oxygen1..mmol.kg., "Water Mass"=metadata$water_mode)
metadata$Sample=parse_number(metadata$Sample)
metadata=metadata[metadata$Sample>=3,]

# Get average values across samples
newmetadata=data.frame(aggregate(list(metadata$Depth, metadata$Cycle, metadata$Temperature, metadata$Salinity, metadata$Oxygen), by=list(metadata$Sample), FUN=mean))
colnames(newmetadata)=c("Sample", "Depth", "Cycle", "Avg Temperature", "Avg Salinity", "Avg Oxygen")

# Relabel water mass types
newmetadata$Water="Subtropical" 
newmetadata$Water[newmetadata$Cycle!=3]="Subantarctic"

# Load coverage data
votus_cov75thres = read.delim(paste0(path, "/data/nz_votus_cov75thres.txt"), header=FALSE, comment.char="#")
colnames(votus_cov75thres) = c("sample","votus","coverage","meandepth","rpkm")

# Load mapping data
viralreads = read.delim(paste0(path,"/data/totalreads_mapped.tsv"), header=FALSE)
totalreads = read.csv(paste0(path, "/data/totalreads.csv"), header=FALSE)
totalreads=totalreads[parse_number(totalreads$V1)>=3,] # Only include samples 3-22 because they have water type labels
viralreads=viralreads[parse_number(viralreads$V1)>=3,]


##### PLOTS #####


### Sequence Length Distribution

# Include only sequences smaller than 50,000 bp - average phage genome length
trimmedseqlens = data.frame("BasePairs" = nzseqlens$V1[nzseqlens$V1<50000])
trimmedseqlens$BasePairs = trimmedseqlens$BasePairs/1000
ggplot(data = trimmedseqlens, aes(x=BasePairs)) +
  geom_histogram(binwidth = .5, color="black", fill="lightblue") +
  xlab("Sequence Length (kBP)") +
  ylab("Number of Sequences") +
  theme_test() +
  ggtitle("Viral Sequence Length Distribution") +
  theme(plot.title = element_text(hjust = .5)) +
  geom_vline(aes(xintercept=mean(BasePairs)), color="black", linetype="dashed", linewidth=.7) +
  annotate("text", x=12, y=1700, label=paste("Mean Sequence \nLength (kBP) = ", round(mean(trimmedseqlens$BasePairs),2)), size=3, hjust = 0)


### Proportion of Viral Reads

propviralreads = data.frame("Sample"=parse_number(viralreads$V1),"Viral"=viralreads$V2/totalreads$V2, "Water"=newmetadata$Water)
propviralreads$Water = factor(propviralreads$Water,levels = c("Subtropical", "Subantarctic"))
colnames(propviralreads) = c("Sample","Viral", "Water")

ggplot(data=propviralreads, aes(y =Viral, x=Water, color=Water)) +
  geom_point() +
  theme_test() +
  xlab("Sample") +
  ylab("Proportion of Reads Mapped") +
  ggtitle("Proportion of Sample Reads Mapped to Viral Sequences")
  


### Total RPKM 

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


# RPKM vs salinity colored by cycle
ggplot(data=tmp,aes(`Avg Salinity`,rpkm,color=Cycle)) +
  #geom_boxplot(color="black") +
  geom_point(size=2.5) +
  theme_test() + 
  xlab("Avg Salinity") + 
  ylab("Total RPKM") +
  theme(axis.text.x=element_text(size=13), axis.text.y=element_text(size=13), legend.text = element_text(size=13), axis.title.x = element_text(size=15),axis.title.y = element_text(size=15))

# RPKM vs temperature colored by cycle
ggplot(data=tmp,aes(`Avg Temperature`,rpkm,color=Cycle)) +
  #geom_boxplot(color="black") +
  geom_point(size=2.5) +
  theme_test() + 
  xlab("Avg Temperature") + 
  ylab("Total RPKM") +
  theme(axis.text.x=element_text(size=13), axis.text.y=element_text(size=13), legend.text = element_text(size=13), axis.title.x = element_text(size=15),axis.title.y = element_text(size=15))


### Number of VOTUS per sample
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


### Shannon diversity for full mapping
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
  ylab("Shannon Diversity index") +
  theme(axis.text.x=element_text(size=13), axis.text.y=element_text(size=13), legend.text = element_text(size=13), axis.title.x = element_text(size=15),axis.title.y = element_text(size=15))



### Cumulative number of vOTUs per sample


perms = 100
le_tmp = table(votus_cov75thres[,1:2])
main = matrix(0,nrow=perms,ncol=20)
for(perm in 1:perms){
  tots = rep(0,18490)
  tmp = le_tmp[sample(1:20),]
  for(i in 1:20){
    tmp2= tmp[i,] - tots
    tmp2[tmp2<0]=0
    main[perm,i] = sum(tmp2)
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


### Rank Order Plot
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
  theme_test() +
  xlab("Order Proportion") +
  ylab("RPKM")


### NMDS Plot
votus_cov75thres=votus_cov75thres[parse_number(votus_cov75thres$sample)>=3,]
otherdf=votus_cov75thres[,c(1,2,5)]
rpkm = dcast(otherdf, sample ~ votus)
rpkm[is.na(rpkm)] = 0
rpkm$sample=parse_number(rpkm$sample)
rpkm=rpkm[rpkm$sample!=19] # Exclude outlier sample
m_rpkm = as.matrix(rpkm)
k5 = kmeans(m_rpkm, 5, nstart=25)
fviz_cluster(k5, m_rpkm, ggtheme = theme_test(), geom = "point")

nmds = metaMDS(m_rpkm, distance="bray")
data.scores = data.frame(scores(nmds)$sites)
data.scores$Sample = rpkm$sample
data.scores=data.scores[order(data.scores$Sample),]
tmp = as.numeric(as.character(votus_cov75thres$rpkm[votus_cov75thres$sample!=19]))
tmp = aggregate(tmp,list(sample=votus_cov75thres$sample[votus_cov75thres$sample!=19]),sum)
colnames(tmp)=c("sample","rpkm")
tmp$sample=parse_number(tmp$sample)
tmp=tmp[tmp$sample>=3,]
tmp=tmp[order(tmp$sample),]
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


### Venn diagram

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

# Can only plot 4 categories at a time, find alternate way of displaying
ggplot(plotcycles, aes(A = as.logical(cycle1), B = as.logical(cycle2), C = as.logical(cycle3), D = as.logical(cycle4), E = as.logical((cycle5)))) +
  geom_venn(show_percentage = FALSE, set_names = c("Cycle 1", "Cycle 2", "Cycle 3", "Cycle 4", "Cycle 5")) + 
  theme_test() +
  theme(text = element_text(size = 15), panel.border = element_blank()) +
  scale_x_continuous(NULL, NULL, NULL) +
  scale_y_continuous(NULL, NULL, NULL)

