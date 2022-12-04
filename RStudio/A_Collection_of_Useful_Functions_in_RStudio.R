### sum(1,2) = 3
sum <- function(a, b){
    return (a + b)
}

has_string_arr <- function(str, arr){
    flag <- 0
    ar <- array()
    for(i in 1:length(arr)){
        if( grepl(str, arr[i], fixed=TRUE) ){
            flag <- flag + 1
            ar[flag] <- i
        }
    }
    return(ar)
}

tri_intersect <- function(A, B, C){
    a <- intersect(A, B)
    b <- intersect(a, C)
    return(b)
}

trim_na <- function(vector){
  for(i in 1: length(vector)){
    if(is.na(vector[i])){
      vector[i] = 0
    }
  }
  return (vector)
}

ignore_na <- function(vector){
  flag = 0
  ar = array()
  for(i in 1: length(vector)){
    if(! is.na(vector[i])){
      flag = flag + 1
      ar[flag] = vector[i]
    }
  }
  return (ar)
}


replace_in_vector <- function(v, from_element, to_element){
  ar <- array(dim = length(v))

  for(i in 1 : length(v)){
    if(v[i] == from_element){
      ar[i] = to_element
    } else{
      ar[i] = v[i]
    }
  }

  return(ar)
}


corner <- function(x, num = 10){
  return(x[1:min(  num, dim(x)[1]  ), 
           1:min(  num, dim(x)[2]  )])
}

my.cbind <- function(data1, data2){
  a <- intersect(rownames(data1), rownames(data2))
  return (cbind(data1[a, ], data2[a, ]))
}

my.rbind <- function(data1, data2){
  a <- intersect(colnames(data1), colnames(data2))
  return (rbind(data1[, a], data2[, a]))
}

eliminate_hyphene <- function(str){
    return(gsub("-", "", str))
}

eliminate_symbol <- function(str){
    return(gsub("(", "", str, fixed = TRUE))
}

eliminate_backward <- function(str){
    return(gsub('x.*$', '', str))
}

eliminate_backward <- function(str){
    return(gsub('\\;.*$', '', str))
}

eliminate_forward <- function(str){
    return(gsub('.*x', '', str))
}

log10_factorial <- function(n){
  if(n == 0){
    return(0)
  }

  out <- 0
  for(i in 1 : n){
    out <- out + log(i) / log(10)
  }
  return(out)
}

my.combination <- function(n, k){
  # return nCk = n! / ((n-k)! k!)
  
  if (n == k || n == 0 || k == 0){
    return(1)
  }

  A = log10_factorial(n)
  B = log10_factorial(n-k)
  C = log10_factorial(k)
  
  log10_nCk = A - B - C
  return(10^(log10_nCk))
}

one_way_2_factor_anova <- function(v1, v2){
  dat <- matrix(0, nrow = ( length(v1) + length(v2) ), ncol = 2 )
  for(i in 1 : length(v1) ){
    dat[i, 1] <- v1[i]
    dat[i, 2] <- 'v1'
  }
  for(i in 1 : length(v2) ){
    dat[i + length(v1), 1] <- v2[i]
    dat[i + length(v1), 2] <- 'v2'
  }
  dat <- as.data.frame(dat)
  
  colnames(dat) <- c('val', 'factor')

  anova_IS <- aov(val ~ factor, data = dat)
  print(summary(anova_IS))

  anova_residuals <- anova_IS$residuals
  print(summary(anova_residuals))
}

one_way_3_factor_anova <- function(v1, v2, v3){
  dat <- matrix(0, nrow = ( length(v1) + length(v2) + length(v3) ), ncol = 2 )
  for(i in 1 : length(v1) ){
    dat[i, 1] <- v1[i]
    dat[i, 2] <- 'v1'
  }
  for(i in 1 : length(v2) ){
    dat[i + length(v1), 1] <- v2[i]
    dat[i + length(v1), 2] <- 'v2'
  }
  for(i in 1 : length(v3) ){
    dat[i + length(v1) + length(v2), 1] <- v3[i]
    dat[i + length(v1) + length(v2), 2] <- 'v3'
  }
  dat <- as.data.frame(dat)
  
  colnames(dat) <- c('val', 'factor')

  anova_IS <- aov(val ~ factor, data = dat)
  print(summary(anova_IS))

  anova_residuals <- anova_IS$residuals
  print(summary(anova_residuals))
}

one_way_4_factor_anova <- function(v1, v2, v3, v4){
  dat <- matrix(0, nrow = ( length(v1) + length(v2) + length(v3) + length(v4) ), ncol = 2 )
  for(i in 1 : length(v1) ){
    dat[i, 1] <- v1[i]
    dat[i, 2] <- 'v1'
  }
  for(i in 1 : length(v2) ){
    dat[i + length(v1), 1] <- v2[i]
    dat[i + length(v1), 2] <- 'v2'
  }
  for(i in 1 : length(v3) ){
    dat[i + length(v1) + length(v2), 1] <- v3[i]
    dat[i + length(v1) + length(v2), 2] <- 'v3'
  }
  for(i in 1 : length(v4) ){
    dat[i + length(v1) + length(v2) + length(v3), 1] <- v4[i]
    dat[i + length(v1) + length(v2) + length(v3), 2] <- 'v4'
  }
  dat <- as.data.frame(dat)
  
  colnames(dat) <- c('val', 'factor')

  anova_IS <- aov(val ~ factor, data = dat)
  print(summary(anova_IS))

  anova_residuals <- anova_IS$residuals
  print(summary(anova_residuals))
}

one_way_5_factor_anova <- function(v1, v2, v3, v4, v5){
  dat <- matrix(0, nrow = ( length(v1) + length(v2) + length(v3) + length(v4) + length(v4) ), ncol = 2 )
  for(i in 1 : length(v1) ){
    dat[i, 1] <- v1[i]
    dat[i, 2] <- 'v1'
  }
  for(i in 1 : length(v2) ){
    dat[i + length(v1), 1] <- v2[i]
    dat[i + length(v1), 2] <- 'v2'
  }
  for(i in 1 : length(v3) ){
    dat[i + length(v1) + length(v2), 1] <- v3[i]
    dat[i + length(v1) + length(v2), 2] <- 'v3'
  }
  for(i in 1 : length(v4) ){
    dat[i + length(v1) + length(v2) + length(v3), 1] <- v4[i]
    dat[i + length(v1) + length(v2) + length(v3), 2] <- 'v4'
  }
  for(i in 1 : length(v5) ){
    dat[i + length(v1) + length(v2) + length(v3) + length(v4), 1] <- v5[i]
    dat[i + length(v1) + length(v2) + length(v3) + length(v4), 2] <- 'v5'
  }
  dat <- as.data.frame(dat)
  
  colnames(dat) <- c('val', 'factor')

  anova_IS <- aov(val ~ factor, data = dat)
  print(summary(anova_IS))

  anova_residuals <- anova_IS$residuals
  print(summary(anova_residuals))
}

comparison_of_two_vectors <- function(v1, v2, paired = FALSE){
  p.val = t.test(v1, v2, paired = paired)
  print(p.val)
  
  log2FC = log( mean(v1 + 0.000000000001)/mean(v2 + 0.000000000001) ) / log(2)
  print(log2FC)
}

my.Fisher.exact.test <- function(total, A, B, cross){
  a1 <- log10_factorial(A)
  a2 <- log10_factorial(total - A)
  a3 <- log10_factorial(B)
  a4 <- log10_factorial(total - B)

  b1 <- log10_factorial(cross)
  b2 <- log10_factorial(A - cross)
  b3 <- log10_factorial(B - cross)
  b4 <- log10_factorial(total - cross - (A - cross) - (B - cross))
  b5 <- log10_factorial(total)

  out = a1 + a2 + a3 + a4 - b1 - b2 - b3 - b4 - b5
  return(10^out)
}

my.MIA.assay.enrichment <- function(total, A, B, cross){
  out <- 0
  for(i in cross:min(A, B)){
    out = out + my.Fisher.exact.test(total, A, B, i)
  }
  return(out)
}

my.MIA.assay.depletion <- function(total, A, B, cross){
  out <- 0
  for(i in 0:cross-1){
    out = out + my.Fisher.exact.test(total, A, B, i)
  }
  return(out)
}

ensembl_to_gene <- function(ensembl_list){
  ar = array(dim = length(ensembl_list))

  human = read.csv("https://blog.kakaocdn.net/dn/BfhRT/btrSyrd9VIx/CDHawSJNTVQs04m292RQs0/human_genes_36601.tsv?attach=1&knm=tfile.tsv", sep = '\t', header = F)
  mouse = read.csv("https://blog.kakaocdn.net/dn/bvsctt/btrSB3cmIyi/AymjWBaekuCgRqWtI4G9Fk/mouse_genes_32285.tsv?attach=1&knm=tfile.tsv", sep = '\t', header = F)

  for(i in 1:length(ensembl_list)){
    if(grepl('ENSG', ensembl_list[i], fixed = TRUE)){ # human gene
      index = match(ensembl_list[i], human[, 1])
      ar[i] = human[index, 2]
    } 
    else if(grepl('ENSMUSG', ensembl_list[i], fixed = TRUE)){ # mouse gene
      index = match(ensembl_list[i], mouse[, 1])
      ar[i] = mouse[index, 2]
    } 
  }
  return(ar)
}

gene_to_ensembl <- function(gene_list){
  ar = array(dim = length(gene_list))

  human = read.csv("https://blog.kakaocdn.net/dn/BfhRT/btrSyrd9VIx/CDHawSJNTVQs04m292RQs0/human_genes_36601.tsv?attach=1&knm=tfile.tsv", sep = '\t', header = F)
  mouse = read.csv("https://blog.kakaocdn.net/dn/bvsctt/btrSB3cmIyi/AymjWBaekuCgRqWtI4G9Fk/mouse_genes_32285.tsv?attach=1&knm=tfile.tsv", sep = '\t', header = F)

  for(i in 1:length(gene_list)){
    if(gene_list[i] == toupper(gene_list[i])){ # human gene
      index = match(gene_list[i], human[, 2])
      ar[i] = human[index, 1]
    } 
    else{ # mouse gene
      index = match(gene_list[i], mouse[, 2])
      ar[i] = mouse[index, 1]
    } 
  }

  # return(ignore_na(ar)) 
  ## if possible, ignore_na should be used

  return(ar)
}

human_to_mouse <- function(human_gene){
  hom <- read.csv("https://blog.kakaocdn.net/dn/OCkAZ/btrL4AnnDXH/YJ0CrvSIchlGTnqvKphemK/HOM_MouseHumanSequence.csv?attach=1&knm=tfile.csv")

  mouse_gene = array()
  flag = 0

  for(i in 1 : length(human_gene)){
    index = match(human_gene[i], hom[hom$Common.Organism.Name == 'human', 'Symbol'])
    key = hom[hom$Common.Organism.Name == 'human', 'DB.Class.Key'][index]
    flag = flag + 1
    mouse_gene[flag] = hom[hom$DB.Class.Key == key 
                           & hom$Common.Organism.Name == 'mouse, laboratory'
                           , 'Symbol'][1] # duplicate mouse genes can be found
  }
	
  return(mouse_gene)
}

mouse_to_human <- function(mouse_gene){
  hom <- read.csv("https://blog.kakaocdn.net/dn/OCkAZ/btrL4AnnDXH/YJ0CrvSIchlGTnqvKphemK/HOM_MouseHumanSequence.csv?attach=1&knm=tfile.csv")

  human_gene = array()
  flag = 0

  for(i in 1 : length(mouse_gene)){
    index = match(mouse_gene[i], hom[hom$Common.Organism.Name == 'mouse, laboratory', 'Symbol'])
    key = hom[hom$Common.Organism.Name == 'mouse, laboratory', 'DB.Class.Key'][index]
    flag = flag + 1
    human_gene[flag] = hom[hom$DB.Class.Key ==  key
                           & hom$Common.Organism.Name == 'human'
                           , 'Symbol'][1] # duplicate human genes can be found
  }

  return(human_gene)
}


# RenameGenesSeurat  ------------------------------------------------------------------------------------
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = HGNC.updated[[i]]$Suggested.Symbol) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA

  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}
# RenameGenesSeurat(obj = SeuratObj, newnames = HGNC.updated.genes)
 
update_cluster_in_seurat_obj <- function(seurat_obj, barcode, cluster){

  # dim(seurat_obj)[1] = length(barcode) = length(cluster)
  mat <- matrix(0, nrow = length(barcode), ncol = 2)
  mat[, 1] = barcode
  mat[, 2] = cluster
  mat = as.data.frame(mat)
  rownames(mat) = barcode

  seurat_obj@meta.data$orig.ident = mat[rownames(seurat_obj@meta.data), 2]
  # you may need to modify the above code
  seurat_obj@active.ident = as.factor(seurat_obj@meta.data$orig.ident)
  
  return (seurat_obj)
}

scatter_plot <- function(x, y, xlab = "x", ylab = "y", point_size = 2, lab_size = 4, png=TRUE){
  library(ggplot2)

  # the lenth(x) must be same with the length(y)
  mat <- matrix(0, nrow = length(x), ncol = 2)
  mat[, 1] = x
  mat[, 2] = y
  colnames(mat) = c(xlab, ylab)
  mat <- as.data.frame(mat)

  if(png){
    png("./scatter_plot.png",width=2000,height=2000,res=500)
    ggplot(mat, aes(x=x, y=y)) + geom_point(shape=19, size=point_size, color="blue") + theme(plot.background = element_blank(),   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(size =1)) +   stat_smooth(method = lm, level=.95, color="grey") + labs(x=xlab, y=ylab, size=lab_size)
    dev.off()
  } else{
    ggplot(mat, aes(x=x, y=y)) + geom_point(shape=19, size=point_size, color="blue") + theme(plot.background = element_blank(),   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(size =1)) +   stat_smooth(method = lm, level=.95, color="grey") + labs(x=xlab, y=ylab, size=lab_size)
  }
}

my.plot <- function(x, y, col){
	# assume that length(x) = length(y) = length(col)

  plot(x_array, y_array, t="n")
  colfunc <- colorRampPalette(c("#000000", "#EB4600", "#FFF800"))
  
  coll = array(dim = length(col))
  for(i in 1 : length(col)){
    coll[i] <- colfunc(100) [as.integer( col[i] / max(col) * 99 + 1)] 
  }
  
  text(x, y, labels = "●", col = coll, cex = 1)
}

conv_spatial_feature_plot <- function(tissue_dir, Tgenes, quality.control = FALSE){
  library(Seurat)
  library(SeuratData)
  library(ggplot2)
  library(cowplot)
  library(dplyr)

  # reference : https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_07_spatial.html
  
  br.sp = Load10X_Spatial(tissue_dir, slice= 'slice1')
  br.sp <- SCTransform(br.sp, assay = "Spatial", verbose = FALSE, variable.features.n = 1000)

  if(quality.control){
    br.sp <- PercentageFeatureSet(br.sp, "^mt-", col.name = "percent_mito")
    br.sp <- PercentageFeatureSet(br.sp, "^Hb.*-", col.name = "percent_hb")
    br.sp <- br.sp[, br.sp$nFeature_Spatial > 500 & br.sp$percent_mito < 25 & br.sp$percent_hb < 20]
  }

  SpatialFeaturePlot(br.sp, features = Tgenes)
}

my.EnhancedVolcano <- function(gene.name, logFC, adj.P.Val, 
                               pCutoff = 0.05, FCcutoff = 0.3,
                               xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5)){
  library(EnhancedVolcano)

  tested <- matrix(0, nrow = length(gene.name), ncol = 2)
  tested <- as.data.frame(tested)
  for(i in 1:length(gene.name)){
    tested[i, 1] <- logFC[i]
    tested[i, 2] <- adj.P.Val[i]
  }
  rownames(tested) <- gene.name
  colnames(tested) <- c('logFC', 'adj.P.Val')

  EnhancedVolcano(tested, lab = rownames(tested), 
                  x='logFC', y='adj.P.Val', xlim = xlim, ylim = ylim, 
                  pCutoff = pCutoff, FCcutoff = FCcutoff) 
}

GO <- function(gene){
    # ont = "ALL", "BP", "CC", "MF"
    # showCategory is not mandatory

    gene <- gsub("GRCh38", "", gene) # human 데이터 가공시의 reference 이름 제거
    gene <- gsub("mm10", "", gene) # mouse 데이터 가공시의 reference 이름 제거
    for(i in 1:10){
  	  gene <- gsub("-", "", gene) # 불필요한 앞부분의 - 제거
    }
    gene <- gsub('\\ .*$', '', gene) # 'KLK2 ENSG00000167751' 같은 것을 해결 
    
    if (gene[1] == toupper(gene[1])){ ## Human gene
        gene.df <- bitr(gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
        gene.df <- as.vector(gene.df[[2]])
        GO <- enrichGO(gene.df, OrgDb = 'org.Hs.eg.db',keyType = "ENTREZID", ont = "ALL", pvalueCutoff = 0.05, pAdjustMethod = "BH")
        return(GO)
	} else{ ## Mouse gene?
        gene.df <- bitr(gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
        gene.df <- as.vector(gene.df[[2]])
        GO <- enrichGO(gene.df, OrgDb = 'org.Mm.eg.db',keyType = "ENTREZID", ont = "ALL", pvalueCutoff = 0.05, pAdjustMethod = "BH")
        return(GO)
    }
}

GO.plot <- function(gene){
  library(EnhancedVolcano)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(enrichplot)

  
    # ont = "ALL", "BP", "CC", "MF"
    # showCategory is not mandatory

    gene <- gsub("GRCh38", "", gene) # human 데이터 가공시의 reference 이름 제거
    gene <- gsub("mm10", "", gene) # mouse 데이터 가공시의 reference 이름 제거
    for(i in 1:10){
  	  gene <- gsub("-", "", gene) # 불필요한 앞부분의 - 제거
    }
    gene <- gsub('\\ .*$', '', gene) # 'KLK2 ENSG00000167751' 같은 것을 해결 
    
    if (gene[1] == toupper(gene[1])){ ## Human gene
        gene.df <- bitr(gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
        gene.df <- as.vector(gene.df[[2]])
        GO <- enrichGO(gene.df, OrgDb = 'org.Hs.eg.db',keyType = "ENTREZID", ont = "ALL", pvalueCutoff = 0.05, pAdjustMethod = "BH")
        dotplot(GO,split="ONTOLOGY", showCategory = 5)+facet_grid(ONTOLOGY~., scale="free")
    } else{ ## Mouse gene?
        gene.df <- bitr(gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
        gene.df <- as.vector(gene.df[[2]])
        GO <- enrichGO(gene.df, OrgDb = 'org.Mm.eg.db',keyType = "ENTREZID", ont = "ALL", pvalueCutoff = 0.05, pAdjustMethod = "BH")
        dotplot(GO,split="ONTOLOGY", showCategory = 5)+facet_grid(ONTOLOGY~., scale="free")
    }
}
