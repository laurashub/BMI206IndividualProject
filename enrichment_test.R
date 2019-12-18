library(ggplot2)
library(xtable)

read_file <- function(filename){
 file_info <- read.table(filename, header = FALSE, sep = " ", dec = ".")
 return(file_info)
}


#Enrichment for essential reactions in top 5% of ranked_data
hypergeom_essential_enrichment <- function(ranked_data){
  essentials <- sum(ranked_data[,3] == 'True')
  nonessentials <- sum(ranked_data[,3] == 'False')

  top_0.05 <- ranked_data[1:floor(dim(ranked_data)[1]*0.05),]
  max_essentials <- dim(top_0.05)[1]
  top_essentials <- sum(top_0.05[,3] == 'True')

  pval <- sum(dhyper(top_essentials:max_essentials, essentials, nonessentials, max_essentials))
  return(pval)
}

#Enrichment for overlapping reactions in top 5% of data1 with data2
hypergeom_overlap_enrichment <- function(data1, data2){ 
  top1_0.05 <- data1[1:floor(dim(data1)[1]*0.05),1]
  top2_0.05 <- data2[1:floor(dim(data2)[1]*0.05),1]
  
  essentials <- sum(sapply(data1, function(x){return(is.element(x, top2_0.05))}))
  nonessentials <- dim(data1)[1] - essentials
  
  max_essentials <- length(top1_0.05)
  top_essentials <- sum(sapply(top1_0.05, function(x){return(is.element(x, top2_0.05))}))
  
  pval <- sum(dhyper(top_essentials:max_essentials, essentials, nonessentials, max_essentials))
  return(pval)
}

#fishers, used to check hypergeometric, returned same values
#fishers_essential_enrichment <- function(ranked_data){
#  top_0.05 <- ranked_data[1:floor(dim(ranked_data)[1]*0.05),]
#  bottom_0.05 <- ranked_data[ceiling(dim(ranked_data)[1]*0.05):dim(ranked_data)[1],]
#  
#  top_essentials <- sum(top_0.05[,3] == 'True')
#  top_nonessentials <- sum(top_0.05[,3] == 'False')
  
#  bottom_essentials <- sum(bottom_0.05[,3] == 'True')
#  bottom_nonessentials <- sum(bottom_0.05[,3] == 'False')
  
#  ftest <- fisher.test(matrix(c(top_essentials,bottom_essentials,top_nonessentials,bottom_nonessentials),
 #             nrow=2,ncol=2))
#  print(ftest)
#}

#read all files into list of tables
strains <- c('MG1655', 'W3110', 'EDL933', 'SAKAI', 'CFTO73', 'UTI89', 'core')
metrics <- c('bridging_centrality', 'betweenness_centrality', 'clustering_coefficient')
data = vector(mode = "list", length = (length(strains) * length(metrics)))
i = 1
for (strain in strains){
  for (metric in metrics){
    name <- sprintf('%s_%s.txt', strain, metric)
    data[[i]] <- read_file(name)
    i <- i+1
  }
}

#calculate essential pvals for each list read in
essential_pvals <- c(1:21)
for (i in 1:21){
  essential_pvals[i] <-hypergeom_essential_enrichment(data[[i]])
}

essential_pvals <- matrix(unlist(essential_pvals), ncol = 3, byrow = TRUE)
essential_pvals <- formatC(essential_pvals, format = "e", digits = 3)
xtable(essential_pvals) #format for latex


overlap_pvals  <- c(1:49) 
index = 1
for (i in 0:6){
  for (j in 0:6){
    if (i == j){
      overlap_pvals[index] = NaN
    } else {
      overlap_pvals[index] <- hypergeom_overlap_enrichment(data[[1+(i*3)]], data[[1+(j*3)]])
    }
    index <- index + 1
  }
}

#make heatmap from overlap_pvals
overlap_pvals.log10 <- unlist(lapply(overlap_pvals, function(x){ -log10(x)}))
overlap_pvals.df <- expand.grid(strain1 = strains, strain2 = strains)
overlap_pvals.df$pval <- overlap_pvals.log10
ggplot(data = overlap_pvals.df, aes(x = strain1, y = strain2)) + geom_tile(aes(fill = pval)) 
ggsave("overlap_pval.png")
