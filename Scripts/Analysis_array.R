
## Salmon array

library(limma)

# Changing to the directory with the files of BlueFuse

setwd("Salmon_raw_data/Liver/")

# the targets file, which contains the files and what is it in each channel.

targets <- readTargets('Targets.txt')

## We will read the files in the target file

RG <- read.maimages(targets, source = 'bluefuse')

# We can see a summary of what we have loaded using show

show(RG)

########## QC control ###########

# We will produce the MAplots
# MAplots should be centered at 0 in the y axis.

plotMA3by2(RG) # this function generates PNGs files in the working directory.

boxplot(data.frame(log2(RG$G)),main="Green background") # the background in the green channel is more or less uniform
                                                        # this is logical, because this is the reference channel


boxplot(data.frame(log2(RG$R)),main="Red background") # The red background is also uniform.




# QC result in this case is bad, we will try to use normalization in order to improve this

########## Normalization #########


# we will try differents methods of normalization in order to improve the data.

# Backgound correction


# Normexp

RG_normexp <- backgroundCorrect(RG, method="normexp", offset=50)


plotMA(RG_normexp)

# # Auto
# 
# RG_auto <- backgroundCorrect(RG, method="auto")
# 
# 
# plotMA(RG_auto)
# 
# # substract
# 
# RG_sub <- backgroundCorrect(RG, method="subtract")
# 
# 
# plotMA(RG_sub)
# 
# # half
# 
# RG_half <- backgroundCorrect(RG, method="half")
# 
# 
# plotMA(RG_half)
# 
# # minimum
# 
# RG_min <- backgroundCorrect(RG, method="minimum")
# 
# 
# plotMA(RG_min)
# 
# # movingmin
# 
# RG_movin <- backgroundCorrect(RG, method="movingmin")
# 
# 
# plotMA(RG_movin)
# 
# # edwards
# 
# RG_edwards <- backgroundCorrect(RG, method="edwards")
# 
# 
# plotMA(RG_edwards)

# We will use the normexp normalization inside each array

# We will use the loess method for ins-array normalization

plotDensities(RG_normexp)

MA <- normalizeWithinArrays(RG_normexp, method="loess")

plotDensities(MA)

MA_norm <- normalizeBetweenArrays(MA, method="Aquantile")

plotDensities(MA_norm) # It's not ideal, but its much better than in the begining

plotMA(MA_norm) # We can see that we have centered the MA plot in 0, solving some problems.


## PCA analysis

# Creation of the design matrix

design <- modelMatrix(targets, ref = 'Ref')

rownames(design) <- c("AC1","AC2","CA1","AT1",
                      "CO1","AC3","CA2","CA3",
                      "AT2","CO2","CA4","AT3",
                      "AT4","CO3","CO4","AC4",
                      "AC5","CA5","AT5","CO5")


samples <- c("AC1","AC2","CA1","AT1",
             "CO1","AC3","CA2","CA3",
             "AT2","CO2","CA4","AT3",
             "AT4","CO3","CO4","AC4",
             "AC5","CA5","AT5","CO5")

colors<-c("red","green","blue","black") 

col_groups <- rep(NA, 20)

col_groups[grep('AC', samples)] <- 'red'

col_groups[grep('CO', samples)] <- 'black'

col_groups[grep('AT', samples)] <- 'blue'

col_groups[grep('CA', samples)] <- 'coral'

library(scatterplot3d)

summary(pca.filt <- prcomp(t(MA_norm$A), scale = TRUE ))

pca3d<-scatterplot3d(x=pca.filt$x[,1],y=pca.filt$x[,2],z=pca.filt$x[,3],
                     xlab='PC1', ylab='PC2', zlab='PC3', main='PCA',
                     pch=16,col.grid="lightblue",
                     color = col_groups)

# text(pca3d$xyz.convert(pca.filt$x+0.5), labels=rownames(design), cex=0.8)## heatmap





################ Differential expression ############


# now we will proceed to the differential expression analysis

fit <- lmFit(MA_norm, design)

# we will define the contrast matrix. We will compare all the compounds against the control group

contrast.matrix <- makeContrasts(AC - Cntrl, AT - Cntrl, CA - Cntrl, levels = design)

fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)

Ac_Cntrl <- topTable(fit2, coef = "AC - Cntrl", number = Inf)

AC_DEGs <- rownames(Ac_Cntrl[Ac_Cntrl$P.Value < 0.01 & abs(Ac_Cntrl$logFC) > 1,])

At_Cntrl <- topTable(fit2, coef = "AT - Cntrl", number = Inf)

At_DEGs <- rownames(At_Cntrl[At_Cntrl$P.Value < 0.01 & abs(At_Cntrl$logFC) > 1,])

Ca_cntrl <- topTable(fit2, coef = "CA - Cntrl", number = Inf)

CA_DEGs <- rownames(Ca_cntrl[Ca_cntrl$P.Value < 0.01 & abs(Ca_cntrl$logFC) > 1,])

################### Venn plots of the results ##########

library(Vennerable)

names <- c('AC', 'AT', 'CA')

list_DEGs <- list(AC_DEGs, At_DEGs, CA_DEGs)

plot(Venn(list_DEGs, SetNames = names))

######### GO annotation ##########

# In this annotation file we have all the DNA sequences of the probes

sequeces <- read.table('cGRASP_16k_sequences.adf.txt', 
                       sep = '\t', header = T,
                       skip = 14, row.names = 1)

# We will obtain the sequences for the AC DEGs and then we will write them to a FASTA file using the seqinr::write.fasta
# function.

AC_seqs <- sequeces[AC_DEGs,]

AC_seqs <- AC_seqs[!(AC_seqs$Reporter.Sequence== ''),] # we wliminate those DEGs that don't have sequence

AC_seqs <- AC_seqs[!(AC_seqs$Reporter.Database.Entry..genbank. == ''),] # and those that dont have gene ID (only one)



seqinr::write.fasta(as.list(AC_seqs$Reporter.Sequence), 
                    names = AC_seqs$Reporter.Database.Entry..genbank.,
                    file.out = 'AC_DEGs.FASTA',open = 'w' )

# We will repeat this porcess for all the DEGs and groups. Now the AT group

AT_seqs <- sequeces[At_DEGs,]

AT_seqs <- AT_seqs[!(AT_seqs$Reporter.Sequence== ''),] # we wliminate those DEGs that don't have sequence

AT_seqs <- AT_seqs[!(AT_seqs$Reporter.Database.Entry..genbank. == ''),] # and those that dont have gene ID 


seqinr::write.fasta(as.list(AT_seqs$Reporter.Sequence), 
                    names = AT_seqs$Reporter.Database.Entry..genbank.,
                    file.out = 'AT_DEGs.FASTA',open = 'w' )

# CA group

CA_seqs <- sequeces[CA_DEGs,]

CA_seqs <- CA_seqs[!(CA_seqs$Reporter.Sequence== ''),] # we wliminate those DEGs that don't have sequence

CA_seqs <- CA_seqs[!(CA_seqs$Reporter.Database.Entry..genbank. == ''),] # and those that dont have gene ID 


seqinr::write.fasta(as.list(CA_seqs$Reporter.Sequence), 
                    names = CA_seqs$Reporter.Database.Entry..genbank.,
                    file.out = 'CA_DEGs.FASTA',open = 'w' )

## Finally we will get the background of the array, all the probes
## so we can use blast2go in order to make an enrichment analysis

background <- sequeces[!(sequeces$Reporter.Sequence == ''),]

background <- background[!(background$Reporter.Database.Entry..genbank. == ''),]

seqinr::write.fasta(as.list(background$Reporter.Sequence),
                    names = background$Reporter.Database.Entry..genbank.,
                    file.out = 'background.fasta', open = 'w')

## Finally we will generate tables containing the candidate genes for
## each group.

library(xlsx)

## CA data

write.xlsx(Ca_cntrl[CA_DEGs,], file = 'CA_candidates.xlsx')

## AT group

write.xlsx(At_Cntrl[At_DEGs,], file = 'AT_candidates.xlsx')

## AC group

write.xlsx(Ac_Cntrl[AC_DEGs,], file = 'AC_candidates.xlsx')

## Heatmap -----

library(gplots)
data.clus <- MA_norm$A
colnames(data.clus) <- samples

clust.rows <- hclust(as.dist(1-cor(t(data.clus))),method="ward")
clust.cols <- hclust(as.dist(1-cor(data.clus)),method="ward.D2")
heatcol<-colorRampPalette(c("green", "Black","red"), space = "rgb")

#pdf(file=file.path(resultsDir, paste("heatmap EC mirna3 unscaled","pdf", sep=".")))

heatm<-heatmap.2(as.matrix(data.clus[1:10000,]), col = heatcol(256),
                 dendrogram="column", Colv=as.dendrogram(clust.cols),
                 Rowv=NULL,
                 scale="row",cexRow=0.1, cexCol=0.5,
                 main="",key=TRUE,keysize=1,density.info="none",trace="none")




