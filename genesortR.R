#Property-based phylogenomic subsampling
#Writen by Nicolas Mongiardino Koch 02/2021

#This script requires an alignment in FASTA format, a partition file in format
#'geneX = 1-200', a species tree considered the best estimate of the true tree
#(e.g., as obtained using concatenation or coalescent methods using the full
#alignment, note that uncertain nodes can be collapsed) and a file with all gene
#trees. Species and gene trees need to be in newick format. The order of genes
#in the alignment must correspond to the order of gene trees order and gene tree
#order need to match. The species tree must be rooted using outgroups.

#Parameters needing input are marked with 'INPUT' and are all in the first
#section below entitled 'Parameters'

#More gene properties than the ones used to infer a usefulness axis are
#inferred. If you would like to sort and subsample based on any of these (such
#as occupancy for example) go to heading 'D) Sort & Subsample and modify the
#lines after 'WARNING'

#More details can be found in the following publications:
#1) Mongiardino Koch & Thompson (2020) - A Total-Evidence Dated Phylogeny of
#Echinoidea Combining Phylogenomic and Paleontological Data. Syst. Biol.
#syaa069, https://doi.org/10.1093/sysbio/syaa069

#2) Mongiardino Koch (2021) - Phylogenomic subsampling and the search for
#phylogenetically reliable loci. bioRxiv 2021.02.13.431075.
#https://doi.org/10.1101/2021.02.13.431075

#Parameters-----------------------------------------------------------------------------------------
#INPUT: set working directory to folder containing these four files
setwd('')

#INPUT: fill this with file names
alignment <- ''
partition <- ''
species_tree <- ''
gene_trees <- ''

#INPUT: is the alignment 'DNA' or 'AA'
type <- 'AA'

#INPUT: provide the names of two terminals that bracket the ingroup, i.e., one
#descendant of each of the two main clades of the ingroup. Leave blank and 
#properties will be claculated across the enitre tree without removing outliers
ingroup <- c('', '')

#INPUT: do not even consider genes with less than 'threshold' ingroup taxa. 
#If threshold == 'auto' then it is automatically set to more than 10% of the
#ingroup terminals (if the dataset is small, a larger value is probably
#desirable)
threshold <- 'auto'

#INPUT: activate/deactivate outlier gene removal (recommended)
remove_outliers <- T
outlier_fraction <- 0.01 #i.e. 1%

##INPUT: Desired number of genes to retain
#if n_genes == 'all' then the dataset is sorted but not subsampled.
n_genes <- 'all'

##INPUT: Whether to incorporate Robinson-Foulds similarity in the PCA. This
##option is available in case the relationships among the studied taxa are
##highly uncertain, and there is concern that specifying a given topology might
##bias results. Alternatively (and preferentially) these uncertainties can also
##be accommodated by using a partially resolved species tree as input, for which
##uncertain relationships have been collapsed. In that case, RF similarities
##will not be affected by the specific resolution of uncertain nodes favored by
##different genes. Note that even if topological_similarity is set to FALSE, a
##species tree needs to be provided to delineate the species in the ingroup
topological_similarity <- T

#Install and load packages-------------------------------------------------------------------------
packages <- c('ape','phytools','phangorn','tibble','dplyr','tidyr','adephylo','ggplot2','cowplot')
new_packages <- packages[!packages %in% installed.packages()[,'Package']]
if(length(new_packages)) { install.packages(new_packages) }

library(phangorn)
library(ape)
library(tibble)
library(dplyr)
library(tidyr)
library(phytools)
library(adephylo)
library(ggplot2)
library(cowplot)

#Some necessary functions-----------------------------------------------------------------------------------
`%not in%` <- function(x, table) is.na(match(x, table, nomatch = NA_integer_))
#function to count invariant sites
inv <- function(x) {
  pattern <- unique(x)
  if(any(pattern %in% c('?'))) pattern <- pattern[-which(pattern == '?')]
  if(any(pattern %in% c('-'))) pattern <- pattern[-which(pattern == '-')]
  if(any(pattern %in% c('X'))) pattern <- pattern[-which(pattern == 'X')]
  
  if(length(pattern) == 1) {
    invariant <- T
  } else {
    invariant <- F
  }
  return(invariant)
}

#function to remove missing data from the estimation of RCFV
remove_empty <- function(x) {
  if('-' %in% names(unlist(x))) {
    missing <- which(names(unlist(x)) == '-')
    x <- x[-missing]
  }
  if('?' %in% names(unlist(x))) {
    missing <- which(names(unlist(x)) == '?')
    x <- x[-missing]
  }
  if('X' %in% names(unlist(x))) {
    missing <- which(names(unlist(x)) == 'X')
    x <- x[-missing]
  }
  return(x)
}

#A) Prepare data-------------------------------------------------------------------------------------------------------
data <- read.phyDat(alignment, format = 'fasta', type = type)
if(type == 'AA') {
  data <- as.AAbin(data) 
} else {
  data <- as.DNAbin(data)
}

partitions <- read.table(partition, sep = ' ')
names <- as.character(unlist(enframe(partitions[,(which(partitions[1,] == '=') - 1)], name = NULL), use.names = F))
partitions <- enframe(partitions[,ncol(partitions)], name = NULL)
partitions <- partitions %>% separate(value, into = c('Start', 'End'), sep = '-') %>% 
  mutate_if(is.character, as.numeric)

gene_trees <- read.tree(gene_trees)
species_tree <- read.tree(species_tree)

#get names of ingroup and outgroup taxa
if(all(nchar(ingroup) != 0)) {
  node <- getMRCA(species_tree, ingroup)
  IG <- Descendants(species_tree, node, type = 'tips')
  IG <- species_tree$tip.label[unlist(IG)]
  OG <- species_tree$tip.label[species_tree$tip.label %not in% IG]
} else {
  IG <- species_tree$tip.label
}

if(threshold == 'auto') {
  threshold <- ceiling(length(IG)/10)
  if(all(nchar(ingroup) != 0)) {
    cat('Setting threshold to evaluate loci to', threshold, 'taxa (i.e., 10% of ingroup taxa).\n')
  } else {
    cat('Setting threshold to evaluate loci to', threshold, 'taxa (i.e., 10% of all taxa).\n')
  }
} else {
  if(is.numeric(threshold)) {
    if(threshold == 0) {
      cat('Taxon threshold is disabled. All loci will be considered regardless of occupancy level.\n')
    } else {
      if(threshold < 1) {
        threshold <- ceiling(length(IG) * threshold)
        cat('Threshold was expecting an integer but was provided a number < 1.\n')
        cat('It will be assumed that this should be taken as a fraction of ingroup taxa.\n')
      } else {
        if(all(nchar(ingroup) != 0)) {
          cat('Loci with less than', threshold, 'ingroup will be discarded.\n')
        } else {
          cat('Loci with less than', threshold, 'taxa will be discarded.\n')
        }
      }
    }
  } else {
    cat("Modify threshold parameter to either \'auto\' or an integer. This run will fail.", '\n')
  }
}

#B) Estimate properties-------------------------------------------------------------------------------
genes <- 1:length(gene_trees)

root_tip_var <- saturation <- missing <- av_patristic <- 
  length <- tree_length <- occupancy <- variable_sites <- 
  RCFV <- rate <- treeness <- average_BS_support <- 
  robinson_sim <- integer(length(gene_trees))

for(i in 1:length(gene_trees)) {
  tree <- gene_trees[[i]]

  #remove genes with less than 'threshold' ingroup taxa
  if(length(which(tree$tip.label %in% IG)) < threshold) next
  
  #if OGs are defined and present in this tree, root with them
  if(length(IG) != length(species_tree$tip.label)) {
    if(any(OG %in% tree$tip.label)) {
      MRCA <- getMRCA(tree, which(tree$tip.label %in% IG))
      tree_rooted <- root(tree, node = MRCA)
      root_tip_var[i] <- var(dist.nodes(tree_rooted)[MRCA, which(tree_rooted$tip.label %in% IG)])
      
      #after this remove the OGs from the gene tree
      tree <- drop.tip(tree, which(tree$tip.label %in% OG))
    } else {
      tree_rooted <- midpoint.root(tree)
      root_tip_var[i] <- var(dist.nodes(tree_rooted)[(length(tree_rooted$tip.label)+1),(1:length(tree_rooted$tip.label))])
    }
  } else {  #otherwise do midpoint rooting
    tree_rooted <- midpoint.root(tree)
    root_tip_var[i] <- var(dist.nodes(tree_rooted)[(length(tree_rooted$tip.label)+1),(1:length(tree_rooted$tip.label))])
  }
  
  average_BS_support[i] <- mean(as.numeric(tree$node.label), na.rm = T)
  if(is.nan(average_BS_support[i])) average_BS_support[i] <- 0
  
  #remove taxa from species tree to match gene tree sampling
  if(length(which(species_tree$tip.label %not in% tree$tip.label)) > 0) {
    this_species_tree <- drop.tip(species_tree, which(species_tree$tip.label %not in% tree$tip.label))
  } else {
    this_species_tree <- species_tree
  }
  
  if(topological_similarity) {
    robinson_sim[i] <- 1 - suppressMessages(RF.dist(this_species_tree, tree, normalize = TRUE, check.labels = TRUE))
  }
  
  patristic_dist <- as.matrix(distTips(tree, tips = 'all', method = 'patristic', useC = T))
  
  #get gene sequence
  gene <- as.character(data[,partitions$Start[i]:partitions$End[i]])
  
  #remove OGs
  if(length(IG) != length(species_tree$tip.label)) {
    gene <- gene[-which(rownames(gene) %in% OG),]
  }
  ntax <- dim(gene)[1]
  
  #remove taxa not in tree (e.g., those with no data for this loci)
  if(any(rownames(gene) %not in% tree$tip.label)) {
    gene <- gene[-which(rownames(gene) %not in% tree$tip.label),]
  }
  
  #remove entirely empty positions (might originate from pruning OGs for
  #example)
  all_missing <- which(apply(gene, 2, function(x) all(x %in% c('-','?','X'))))
  if(length(all_missing) > 0) {
    gene <- gene[,-all_missing]
  }
  
  variable_sites[i] <- 1 - (length(which(apply(gene, 2, inv)))/dim(gene)[2])
  missing[i] <- length(which(gene %in% c('-','?','X')))/(dim(gene)[1]*dim(gene)[2])
  length[i] <- dim(gene)[2]
  occupancy[i] <- dim(gene)[1]/ntax
  
  p_dist <- as.matrix(dist.hamming(as.phyDat(gene, type = type), ratio = TRUE, exclude = "pairwise"))
  p_dist <- p_dist[order(colnames(p_dist)),order(colnames(p_dist))]
  patristic_dist <- patristic_dist[order(colnames(patristic_dist)),order(colnames(patristic_dist))]
  p_dist <- p_dist[lower.tri(p_dist)]
  patristic_dist <- patristic_dist[lower.tri(patristic_dist)]
  av_patristic[i] <- mean(patristic_dist)
  saturation[i] <- 1 - lm(p_dist ~ patristic_dist)$coefficients[[2]]
  if(is.na(saturation[i])) saturation[i] <- 0
  if(saturation[i] > 1) saturation[i] <- 1
  if(saturation[i] < 0) saturation[i] <- 0
  
  #if sequence is made of AA, calculate comp. heterogeneity
  if(type == 'AA') {
    mean_freqs <- table(c(gene))
    states <- sort(unique(c(gene)))
    if('-' %in% states) {
      mean_freqs <- mean_freqs[-which(states == '-')]
      states <- states[-which(states == '-')]
    }
    if('?' %in% states) {
      mean_freqs <- mean_freqs[-which(states == '?')]
      states <- states[-which(states == '?')]
    }
    if('X' %in% states) {
      mean_freqs <- mean_freqs[-which(states == 'X')]
      states <- states[-which(states == 'X')]
    }
    mean_freqs <- mean_freqs/sum(mean_freqs)
    
    freqs <- lapply(split(gene, seq(nrow(gene))), table)
    freqs <- lapply(freqs, remove_empty)
    for(j in 1:length(freqs)) {
      if(!all(states %in% names(freqs[[j]]))) {
        miss <- states[which(states %not in% names(freqs[[j]]))]
        add <- rep(0, length(miss))
        names(add) <- miss
        add <- as.table(add)
        freqs[[j]] <- as.table(c(freqs[[j]], add))
        freqs[[j]] <- freqs[[j]][order(factor(names(freqs[[j]])))]
      }
    }
    
    freqs <- lapply(freqs, function(x) x/sum(x))
    freqs <- lapply(freqs, function(x) abs(x-mean_freqs))
    freqs <- lapply(freqs, function(x) x/dim(gene)[1])
    freqs <- lapply(freqs, function(x) sum(unlist(x)))
    RCFV[i] <- sum(unlist(freqs))
  }
  
  tree_length[i] <- sum(tree$edge.length)
  rate[i] <- sum(tree$edge.length)/length(tree$tip.label)
  treeness[i] <- 1 - (sum(tree$edge.length[which(tree$edge[,2] %in% c(1:length(tree$tip.label)))])/sum(tree$edge.length))
  if(is.nan(treeness[i])) treeness[i] <- 0
}

#gather gene properties
variables <- data.frame(genes, root_tip_var, saturation, missing, rate, tree_length, treeness, 
                        av_patristic, RCFV, length, occupancy, variable_sites, average_BS_support, robinson_sim)



if(type == 'DNA') {
  variables <- variables[,-which(colnames(variables) == 'RCFV')]
} else {
  if(any(is.na(variables$RCFV))) {
    cat('This script will soon crash, as loci and gene trees are not composed of the same taxa.\n')
    cat('Most likely this is due to loci and gene trees not being in the same order.\n')
    cat('This makes the estimation of some gene properties to fail\n')
  }
}

if(!topological_similarity) {
  variables <- variables[,-which(colnames(variables) == 'robinson_sim')]
}
  
#remove those with less than 'threshold' taxa
useless <- which(apply(variables[,-1], 1, function(x) all(x == 0)))
if(length(useless) > 0) {
  variables <- variables[-useless,]
}

#Select gene properties for PCA (RCFV is only included if type == 'AA)
variables_to_use <- which(colnames(variables) %in% c('root_tip_var', 'saturation', 'av_patristic', 'RCFV', 
                                                     'variable_sites', 'average_BS_support', 'robinson_sim'))

if(any(is.na(variables[,variables_to_use]))) {
  cat('Something went wrong. Most likely the order of genes and gene trees does not match.\n')
}

#perform PCA
PCA <- princomp(variables[,variables_to_use], cor = T, scores = T)

column_names <- c('Gene name', 'Position in dataset', 'Root-tip var.', 
                  'Saturation', 'Missing data', 'Evolutionary rate', 
                  'Tree length', 'Treeness', 'Av patristic dist.',
                  'Comp. heterogeneity', 'Alignment length', 'Occupancy', 
                  'Prop. variable sites', 'Av. boostrap support', 
                  'RF similarity')

if(type == 'DNA') column_names <- column_names[-10]

if(remove_outliers) {
  #estimate Mahalanobis distances
  maha_distances <- order(mahalanobis(PCA$scores, rep(0, length(variables_to_use)), cov(PCA$scores)), decreasing = T)
  if(outlier_fraction >= 1/nrow(variables)) {
    #remove the loci within the top outlier_fraction
    outliers <- maha_distances[1:floor(nrow(variables)*outlier_fraction)]
    outliers <- sort(outliers)
    
    outlier_properties <- data.frame(variables[outliers,])
    outlier_properties <- cbind(data.frame(names = names[outliers]), outlier_properties)
    
    colnames(outlier_properties) <- column_names[1:ncol(outlier_properties)]
    write.csv(outlier_properties, file = paste0(getwd(), '/properties_outliers.csv'), row.names = F)
    
    #redo PCA
    PCA <- princomp(variables[-outliers,variables_to_use], cor = T, scores = T)
    
    #get scores for loci along dimensions 1 and 2
    PC_1 <- PCA$scores[,1]
    PC_2 <- PCA$scores[,2]
    for(i in 1:length(outliers)) {
      if(outliers[i] == 1) {
        PC_1 <- c(NA, PC_1)
        PC_2 <- c(NA, PC_2)
      } else {
        if(outliers[i] < length(PC_1)) {
          PC_1 <- c(PC_1[1:(outliers[i]-1)], NA, PC_1[outliers[i]:length(PC_1)])
          PC_2 <- c(PC_2[1:(outliers[i]-1)], NA, PC_2[outliers[i]:length(PC_2)])
        } else {
          PC_1 <- c(PC_1, NA)
          PC_2 <- c(PC_2, NA)
        }
      }
    }
  } else {
    cat('Not enough genes to remove even 1 loci, use different outlier_fraction.\n')
    PC_1 <- PCA$scores[,1]
    PC_2 <- PCA$scores[,2]
  }
} else {
  PC_1 <- PCA$scores[,1]
  PC_2 <- PCA$scores[,2]
}

variables <- cbind(variables, PC_1, PC_2)
if(any(is.na(variables))) variables <- variables[which(complete.cases(variables)),]

#Deal with subsampling threshold
if(n_genes == 'all') {
  n_genes <- nrow(variables)
  cut <- F
} else {
  if(is.character(n_genes)){
    n_genes <- as.numeric(n_genes)
  }
  cut <- T
}

#C) Attempt to find usefulness axis automatically----------------------------------------------------

#if correlation between PC1 and rate is high
if(cor.test(variables$rate, variables$PC_1)$estimate > 0.7) { 
  #and much higher (x2) than correlation between PC2 and rate
  if(cor.test(variables$rate, variables$PC_1)$estimate > (cor.test(variables$rate, variables$PC_2)$estimate * 2)) {
    PC_rate <- 'PC_1'
  } else {
    PC_rate <- 'PC_1_maybe'
  }
} else { #try to see if rate is captured by PC2
  if(cor.test(variables$rate, variables$PC_2)$estimate > 0.7) { 
    #and much higher (x2) than correlation between PC2 and rate
    if(cor.test(variables$rate, variables$PC_2)$estimate > (cor.test(variables$rate, variables$PC_1)$estimate * 2)) {
      PC_rate <- 'PC_2'
    } else {
      PC_rate <- 'PC_2_maybe'
    }
  } else {
    PC_rate <- 'unknown'
  }
}

if(PC_rate != 'unknown') {
  #is rate also usefulness?? i.e. should we choose the fastest evolving loci?
  loadings_usefulness <- loadings(PCA)[][,as.numeric(unlist(strsplit(PC_rate, '_'))[2])]
  biases <- which(names(loadings_usefulness) %in% c('root_tip_var', 'saturation', 'av_patristic', 'RCFV'))
  signal <- which(names(loadings_usefulness) %in% c('average_BS_support', 'robinson_sim'))
  if(all(loadings_usefulness[signal] < 0)) {
    if(length(which(loadings_usefulness[biases] > 0)) >= 2) {
      PC_usefulness <- as.numeric(unlist(strsplit(PC_rate, '_'))[2])
      direction <- 'clear'
      descending <- F
    } else {
      direction <- 'unclear'
    }
  } else {
    if(all(loadings_usefulness[signal] > 0)) {
      if(length(which(loadings_usefulness[biases] < 0)) >= 2) {
        PC_usefulness <- as.numeric(unlist(strsplit(PC_rate, '_'))[2])
        direction <- 'clear'
        descending <- T
      } else {
        direction <- 'unclear'
      }
    } else {
      direction <- 'unclear'
    }
  }
  
  #rate != usefulness, can we find usefulness??
  if(direction == 'unclear') {
    if(as.numeric(unlist(strsplit(PC_rate, '_'))[2]) == 1) {
      PC_usefulness <- 2
    } else {
      PC_usefulness <- 1
    }
    
    loadings_usefulness <- loadings(PCA)[][,PC_usefulness]
    if(all(loadings_usefulness[signal] < 0)) {
      if(length(which(loadings_usefulness[biases] > 0)) >= 2) {
        direction <- 'clear'
        descending <- F
        cat('A usefulness axis has been found!\n')
      } else {
        direction <- 'unclear'
      }
    } else {
      if(all(loadings_usefulness[signal] > 0)) {
        if(length(which(loadings_usefulness[biases] < 0)) >= 2) {
          direction <- 'clear'
          descending <- T
          cat('A usefulness axis has been found!\n')
        }
        direction <- 'unclear'
      } else {
        direction <- 'unclear'
      }
    }
  } else {
    cat('Rate == usefulness. Proceed and loci will be sorted by rate.\n')
  }
}

#D) Sort & Subsample------------------------------------------------------------------------------
if(grepl('maybe', PC_rate)) {
  cat('There seems to be some ambiguity as to the identity of the axes.\n') 
  cat('Proceed with caution and check PCA loadings to see if sorting is appropriate.\n')
}

if(direction == 'unclear') {
  cat('It is unclear how to sort the data.\n') 
  cat('You can the check loadings and decide manually how to proceed.\n') 
  cat('In the absense of a clear usefulness axis my best guess is to sort by rates.\n')
  cat('(i.e., choose the slowest evolving genes)\n')
  
  variables_sorted <- variables[order(variables[,'rate'], decreasing = F),]
}

if(direction == 'clear') {
  cat('Usefulness axis explains', round((PCA$sdev[PC_usefulness]^2/sum(PCA$sdev^2))*100, digits = 2), 
      'percentage of variance\n')
  
  usefulness_col <- grep(PC_usefulness, colnames(variables))
  
  #sort by usefulness
  variables_sorted <- variables[order(variables[,usefulness_col], decreasing = descending),]
}

#output properties of the dataset
variables_sorted_tosave <- cbind(data.frame(names = names[variables_sorted$genes]), variables_sorted)
if(topological_similarity) {
  colnames(variables_sorted_tosave) <- c(column_names, 'PC1', 'PC2')
} else {
  colnames(variables_sorted_tosave) <- c(column_names[-length(column_names)], 'PC1', 'PC2')
}

write.csv(variables_sorted_tosave, file = paste0(getwd(), '/properties_sorted_dataset.csv'), row.names = F)

###WARNING: uncomment and omdify the following lines if you would like to sort
###and subsample by a different property, for example occupancy or RF similarity
#variables_sorted <- variables[order(variables[,'occupancy'], decreasing = T),]
#variables_sorted <- variables[order(variables[,'robinson_sim'], decreasing = T),]

#sort entire dataset according to the sorting order imposed
positions <- c()
for(j in 1:nrow(variables_sorted)) {
  positions <- c(positions, partitions$Start[variables_sorted$genes[j]]:partitions$End[variables_sorted$genes[j]])
}
sorted_data <- data[,positions]

#sort partitions
sorted_partitions <- partitions[variables_sorted$genes,]
for(j in 1:nrow(sorted_partitions)) {
  dif <- sorted_partitions$End[j]-sorted_partitions$Start[j]
  if(j == 1) {
    sorted_partitions$Start[j] <- 1
  } else {
    sorted_partitions$Start[j] <- sorted_partitions$End[j-1]+1
  }
  sorted_partitions$End[j] <- sorted_partitions$Start[j]+dif
}

#sort gene names
sorted_names <- names[variables_sorted$genes]

#sort genes
sorted_trees <- gene_trees[variables_sorted$genes]

#subsample (or not if n_genes was left == 0)
sorted_names <- sorted_names[1:n_genes]
sorted_partitions <- sorted_partitions[1:n_genes,]
sorted_data <- sorted_data[,1:sorted_partitions$End[n_genes]]
sorted_trees <- sorted_trees[1:n_genes]

write.phyDat(phyDat(sorted_data, type = type), 
             file = paste0(getwd(), '/sorted_alignment_', 
                           n_genes, 'genes.fa'), format = 'fasta')
partitions_tosave <- paste0(sorted_names, ' = ', sorted_partitions$Start, '-', sorted_partitions$End)
write(partitions_tosave, file = paste0(getwd(), '/sorted_alignment_', n_genes, 'genes.txt'))
write.tree(sorted_trees, file = paste0(getwd(), '/sorted_trees_', n_genes, 'genes.tre'))

#E) Optional: visualize some sorting results------------------------------------------------------------------------------------
variables_to_plot <- data.frame(gene = rep(variables_sorted$genes, ncol(PCA$loadings)),
                               value = c(as.matrix(variables_sorted[,variables_to_use])),
                               property = rep(colnames(variables_sorted[,variables_to_use]), 
                                              each = nrow(variables_sorted)),
                               pos = rep(1:nrow(variables_sorted), ncol(PCA$loadings)))

order_properties <- as.character(unique(variables_to_plot$property))

if(type == 'DNA') {
  if(topological_similarity) {
    labs <- c('Root-to-tip variance', 'Level of saturation', 'Av. patristic distance', 
              'Prop. of variable sites', 'Average bootstrap', 'RF similarity')
    colors <- c('#8A2B0E', '#C75E24', '#C69E57', '#868568', '#5F7881', '#586160')
  } else {
    labs <- c('Root-to-tip variance', 'Level of saturation', 'Av. patristic distance', 
              'Prop. of variable sites', 'Average bootstrap')
    colors <- c('#8A2B0E', '#C75E24', '#C69E57', '#868568', '#5F7881')
  }
} else {
  if(topological_similarity) {
    labs <- c('Root-to-tip variance', 'Level of saturation', 'Av. patristic distance', 'Comp. heterogeneity', 
              'Prop. of variable sites', 'Average bootstrap', 'RF similarity')
    colors <- c('#5F1202', '#8A2B0E', '#C75E24', '#C69E57', '#868568', '#5F7881', '#586160')
  } else {
    labs <- c('Root-to-tip variance', 'Level of saturation', 'Av. patristic distance', 'Comp. heterogeneity', 
              'Prop. of variable sites', 'Average bootstrap')
    colors <- c('#5F1202', '#8A2B0E', '#C75E24', '#C69E57', '#868568', '#5F7881')
  }
}

for(i in 1:length(unique(variables_to_plot$property))) {
  main_plot <- ggplot(subset(variables_to_plot, property == as.character(unique(variables_to_plot$property)[i])), 
                      aes(x = pos, y = value, color = property)) + 
    geom_point(alpha = 0.2, shape = 16) + geom_smooth(method = 'gam', se = F) +
    theme_bw() + theme(legend.position = "none") +
    xlab('Sorted position') + ylab('Value') + scale_color_manual(values = colors[i]) + ggtitle(labs[i]) + 
  theme(plot.title = element_text(hjust = 0.5))
  
  inset_plot <- ggplot(subset(variables_to_plot, property == as.character(unique(variables_to_plot$property)[i])), 
                       aes(x = pos, y = value, color = property)) + geom_smooth(method = 'gam', se = T) +
    theme_bw() + theme(legend.position = "none") + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
          axis.text.x = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), axis.ticks.x = element_blank()) + 
    scale_color_manual(values = colors[i])
  
  if(cut) inset_plot <- inset_plot + geom_vline(xintercept = n_genes, linetype = 'dashed')
  
  if(length(unique(variables_to_plot$property)) == 7) second_row = 5
  if(length(unique(variables_to_plot$property)) <= 6) second_row = 4
  
  if(i < second_row) {
    plot_with_inset <- ggdraw() + suppressMessages(draw_plot(main_plot)) + 
      suppressMessages(draw_plot(inset_plot, x = 0.15, y = 0.67, width = .3, height = .25))
  } else {
    plot_with_inset <- ggdraw() + suppressMessages(draw_plot(main_plot)) + 
      suppressMessages(draw_plot(inset_plot, x = 0.15, y = 0.10, width = .3, height = .25))
  }
  
  assign(paste0('plot', letters[i]), plot_with_inset)
}

if(type == 'AA') {
  if(topological_similarity) {
    final_plot <- plot_grid(plota, plotb, plotc, plotd, 
                            plote, plotf, plotg, nrow = 2)
  } else {
    final_plot <- plot_grid(plota, plotb, plotc, plotd, 
                            plote, plotf, nrow = 2)
  }
} else {
  if(topological_similarity) {
    final_plot <- plot_grid(plota, plotb, plotc, plotd, 
                            plote, plotf, nrow = 2)
  } else {
    final_plot <- plot_grid(plota, plotb, plotc, plotd, 
                            plote, nrow = 2)
  }
}

#plot and save
plot(final_plot)
ggsave(paste0('sorted_figure_', n_genes, 'genes.pdf'), plot = final_plot, 
       width = 16, height = 9, units = 'in')