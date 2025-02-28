
#####################################################################
#  
#  If working on Great Lakes, be sure you requested "load Bioinformatics Rchipenrich" in the Modules section!

#  If working on your laptop, first install the chipenrich Bioconductor package with the following code.
# install.packages("BiocManager")  # You may not need this
# BiocManager::install("chipenrich")  # this may take several minutes.

library(chipenrich)  # Load the chipenrich package. May take ~20-30 seconds.

?chipenrich  # See a description of the original GSE function. This works well for most transcription factor peak sets.
# This function models whether or not each gene has at least one genome region assigned to it.

?polyenrich  # This is an alternative function that works better for more than ~60,000 genomic regions.
# This function models the number of genomic regions each gene has assigned to it.

supported_genomes()  # See which reference genomes are supported

supported_locusdefs() # See the choices for assigning peaks to genes. 
?supported_locusdefs  # See descriptions of the choices

supported_genesets() # See the choices for annotation databases
# GOBP is Gene Ontology biological processes, etc.

supported_methods() # See that the package contains the functions:
# broadenrich(), chipenrich(), polyenrich()
# For Fisher's exact test (fet), you would use chipenrich(...,method="fet")

# If you're working on your own laptop, for convenience set your working directory to where you saved the peak file.
#  setwd(" --Your path here-- /Lab-GSE")
getwd()


# We will now test what Panther pathways are regulated by JUNB using
# polyenrich, which models the number of peaks per gene.

# First, we will test regulation of pathways from across the genome (nearest TSS), and then
# we will test regulation from promoter regions only (<5kb from a TSS). 
## Normally, one would probably want to test all of GO and some other pathways,
##   but we're just testing Panther pathways now to finish quickly for lab.

# The path used below is for Great Lakes. You will need to change it if working on your laptop
pe.nearestTSS<-polyenrich("/nfs/turbo/dcmb-class/bioinf545/shared/lab4/JunB-64781_peaks.bed",genesets="panther_pathway",
                       genome="hg38", locusdef="nearest_tss")

# setwd(c("C:\\Users/sartorma/Dropbox (University of Michigan)/from_box/Documents/teaching/Bioinf545/W24/Lab-GSE/"))

# Note that you will get some warnings if there are strange chromosomes in the list, which there are in this case. You can ignore these.

# Notice that the following results were saved in your working directory:
#  1.) A *_qcplots.pdf file containing QC plots, including distribution of peaks around TSSs
#  2.) A *_results file containing the enrichment results
#  3.) A *_peaks file showing which gene each peak is assigned to and its distance
#  4.) A *_peaks-per-gene file containing the number of peaks for each gene and each gene's locus length.
#  5.) A *_opts file containing the options used in the test

attributes(pe.nearestTSS)  # These same items are in the R object of results
pe.nearestTSS$results[1:10,]  # See top 10 significant results
# How many pathways are significantly enriched at the FDR < 0.05 level?

# Notice the "status" column indicates whether each pathway was enriched or depleted. 
#  Were the top pathways enriched or depleted?  
#  A significant depleted pathway would mean there are fewer peaks there than expected by chance.

# Next we'll test enrichment using only the peaks within 5kb of a TSS
pe.5kb<-polyenrich("/nfs/turbo/dcmb-class/bioinf545/shared/lab4/JunB-64781_peaks.bed",out_name="polyenrich-5kb",genesets="panther_pathway",
      genome="hg38", locusdef="5kb")

pe.5kb$results[1:5,]

##  How do the results compare?  
##  Is the same pathway most significant in both analyses?

##  Were there any pathways significant for one analysis, but not the other?
##  How would you interpret a pathway that was significant using nearest TSS, but not for 5kb?
##  Merge them and plot correlation of significance.

pe.both= merge(x=pe.nearestTSS$results,y=pe.5kb$results,
               by=c("Geneset.Type","Geneset.ID","Description"))
pe.both[1:2,]
plot(-log10(pe.both$P.value.x),-log10(pe.both$P.value.y))
lines(c(0,6),c(0,6),col="red")
# The above plot shows the significance of enrichment for all peaks using the nearest TSS (x-axis)
#  versus only peaks within 5kb of a transcription start site (y-axis)


############ Create Dot Plot using ggplot2  ##################
##  This code create a dotplot of the top 10 enriched and depleted for both nearest TSS and <5kb
library(ggplot2)
#library(org.Hs.eg.db)

head(pe.nearestTSS$results[,1:11])
head(pe.5kb$results[,1:11])
nearestTSS<-pe.nearestTSS$results[,c(3,4,5,8,9,10)]
results5kb<- pe.5kb$results[,c(3,4,5,8,9,10)]

num_pathway=10
### Test 1: nearest TSS
up1<-nearestTSS[nearestTSS$Status=="enriched",]
up1<-up1[1:num_pathway,] 
down1<-nearestTSS[nearestTSS$Status=="depleted",]
down1<-down1[1:num_pathway,]
all1 <- rbind.data.frame(up1, down1)
all1$method <- c(rep("Enriched", num_pathway), rep("Depleted", num_pathway))
all1

### Test 2: Control Fast Agers vs. Control Slow Agers
up2<-results5kb[results5kb$Status=="enriched",]
up2<-up2[1:num_pathway,] 
down2<-results5kb[results5kb$Status=="depleted",]
down2<-down2[1:num_pathway,]
all2 <- rbind.data.frame(up2, down2)
all2$method <- c(rep("Enriched", num_pathway), rep("Depleted", num_pathway))
all2

### Combine Tests
all <- rbind.data.frame(all1, all2)
all$comparison <- c(rep("nearest TSS", 2*num_pathway),
                    rep("5kb", 2*num_pathway))
all$comparison<-factor(all$comparison, levels=c("nearest TSS","5kb"))
#colnames(all)[7] <- "FDR"
all$Description<-factor(all$Description,levels=unique(all$Description))

######## Dotplot
#png("dotplot_GSEA.png", 
 #   height = 9, width = 7, units = "in", res = 300)
#  If you want to save to a file instead of outputting to the screen, you can run the png() function and the dev.off() function at the end
ggplot(all, 
       aes(x=factor(method, levels = c("Enriched","Depleted")),
           y=Description)) + 
  geom_point(aes(size=abs(N.Geneset.Genes), color=FDR)) +
  scale_colour_gradient(low="red", high="darkorchid1") +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle=45,hjust=0.9,size=12),
        panel.background = element_rect(fill = "white",
                                        colour = "white"),
        panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    linewidth = 1),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "lightgrey")) +
  scale_y_discrete(limits=rev) +
  facet_wrap(~comparison)
#dev.off()

# Finally, you may want to write the results to a file, do additional filtering, create visualizations, etc.

# !!!If using Great Lakes, When finished be sure to log out and close your job so you don't waste any credits. 
