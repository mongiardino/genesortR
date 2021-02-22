#Property-based phylogenomic subsampling
#Writen by Nicolas Mongiardino Koch 02/2021

#This script requires an alignment in FASTA format, a partition file in format
#'geneX = 1-200', a species tree considered the best estimate of the true tree
#(e.g., as obtained using concatenation or coalescent methods using the full
#alignment, not that uncertain nodes can be collapsed) and file with all gene
#trees. Species and gene trees need to be in newick format. The order of genes
#in the alignment must correspond to the order of gene trees order and gene tree
#order need to match. The species tree must be rooted using outgroups.

#Parameters needing input are marked with 'INPUT' and are all in the second
#section called 'Parameters'

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

#Install and load packages-------------------------------------------------------------------------
packages <- c('ape','phytools','phangorn','tibble','dplyr','tidyr','adephylo','ggplot2')
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

#Parameters-----------------------------------------------------------------------------------------
genesortR <- function()
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
#decendant of each of the two main clades of the ingroup
ingroup <- c('', '')

#INPUT: do not even consider genes with less than 'threshold' ingroup taxa
threshold <- 5

#INPUT: activate/deactivate outlier gene removal (recommended)
remove_outliers <- T
outlier_fraction <- 0.01 #i.e. 1%

##INPUT: Desired number of genes to retain
#if n_genes == 0 (i.e. if next line is not modified, then the dataset is sorted
#and saved without subsampling)
n_genes <- 0

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

if(n_genes == 0) {
  n_genes <- nrow(partitions)
}

gene_trees <- read.tree(gene_trees)
species_tree <- read.tree(species_tree)

#get names of ingroup and outgroup taxa
node <- getMRCA(species_tree, ingroup)
IG <- Descendants(species_tree, node, type = 'tips')
IG <- species_tree$tip.label[unlist(IG)]
OG <- species_tree$tip.label[species_tree$tip.label %not in% IG]

#B) Estimate properties-------------------------------------------------------------------------------
root_tip_var <- saturation <- missing <- av_patristic <- length <- tree_length <- occupancy <- variable_sites <- RCFV <- 
  rate <- average_BS_support <- robinson_sim <- vector(length = length(gene_trees))

for(i in 1:length(gene_trees)) {
  tree <- gene_trees[[i]]
  #remove genes with less than 'threshold' ingroup taxa
  if(length(which(tree$tip.label %in% IG)) < threshold) next
  
  #if OGs are present, root with them
  if(any(OG %in% tree$tip.label)) {
    MRCA <- getMRCA(tree, which(tree$tip.label %in% IG))
    tree_rooted <- root(tree, node = MRCA)
    root_tip_var[i] <- var(dist.nodes(tree_rooted)[MRCA, which(tree_rooted$tip.label %in% IG)])
    
    #otherwise do midpoint rooting
  } else {
    tree_rooted <- midpoint.root(tree)
    root_tip_var[i] <- var(dist.nodes(tree_rooted)[(length(tree_rooted$tip.label)+1),(1:length(tree_rooted$tip.label))])
  }
  
  tree <- drop.tip(tree, which(tree$tip.label %in% OG))
  average_BS_support[i] <- mean(as.numeric(tree$node.label), na.rm = T)
  
  this_species_tree <- drop.tip(species_tree, which(species_tree$tip.label %not in% tree$tip.label))
  robinson_sim[i] <- 1 - suppressMessages(RF.dist(this_species_tree, tree, normalize = TRUE, check.labels = TRUE))
  patristic_dist <- as.matrix(distTips(tree, tips = 'all', method = 'patristic', useC = T))
  
  #get gene sequence
  gene <- as.character(data[,partitions$Start[i]:partitions$End[i]])
  
  #remove OGs
  gene <- gene[-which(rownames(gene) %in% OG),]
  ntax <- dim(gene)[1]
  
  #remove taxa not in tree (e.g., those with no data for this loci)
  if(any(rownames(gene) %not in% tree$tip.label)) {
    gene <- gene[-which(rownames(gene) %not in% tree$tip.label),]
  }
  
  #remove entirely empty positions (might originate from prunnin OGs for
  #example)
  all_missing <- which(apply(gene, 2, function(x) all(x %in% c('-','?'))))
  if(length(all_missing) > 0) {
    gene <- gene[,-all_missing]
  }
  
  variable_sites[i] <- 1 - length(which(apply(gene, 2, inv)))/dim(gene)[2]
  missing[i] <- length(which(gene %in% c('-','?')))/(dim(gene)[1]*dim(gene)[2])
  length[i] <- dim(gene)[2]
  occupancy[i] <- dim(gene)[1]/ntax
  
  p_dist <- as.matrix(dist.hamming(as.phyDat(gene, type = type), ratio = TRUE, exclude = "pairwise"))
  p_dist <- p_dist[order(colnames(p_dist)),order(colnames(p_dist))]
  patristic_dist <- patristic_dist[order(colnames(patristic_dist)),order(colnames(patristic_dist))]
  p_dist <- p_dist[lower.tri(p_dist)]
  patristic_dist <- patristic_dist[lower.tri(patristic_dist)]
  av_patristic[i] <- mean(patristic_dist)
  saturation[i] <- 1 - lm(p_dist ~ patristic_dist)$coefficients[[2]]
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
}

#gather gene properties
variables <- data.frame(genes = 1:length(gene_trees), root_tip_var, saturation, missing, rate, tree_length, 
                        RCFV, av_patristic, length, occupancy, variable_sites, average_BS_support, robinson_sim)
if(type == 'DNA') variables <- variables[,-which(colnames(variables) == 'RCFV')]

#remove those with less than 'threshold' taxa
useless <- which(apply(variables[,-1], 1, function(x) all(x == 0)))
if(length(useless) > 0) {
  variables <- variables[-useless,]
}

#Select gene properties for PCA (RCFV is only included if type == 'AA)
variables_to_use <- which(colnames(variables) %in% c('root_tip_var', 'saturation', 'RCFV', 'av_patristic', 
                                                    'variable_sites', 'average_BS_support', 'robinson_sim'))
#perform PCA
PCA <- princomp(variables[,variables_to_use], cor = T, scores = T)

if(remove_outliers) {
  #estimate Mahalanobis distances
  maha_distances <- order(mahalanobis(PCA$scores, rep(0, length(variables_to_use)), cov(PCA$scores)), decreasing = T)
  if(outlier_fraction >= 1/nrow(variables)) {
    #remove the loci within the top outlier_fraction
    outliers <- maha_distances[1:floor(nrow(variables)*outlier_fraction)]
    outliers <- sort(outliers)
    
    #redo PCA
    PCA <- princomp(variables[-outliers,variables_to_use], cor = T, scores = T)
    
    #get scores for loci along dimensions 1 and 2
    scores1 <- PCA$scores[,1]
    scores2 <- PCA$scores[,2]
    for(i in 1:length(outliers)) {
      if(outliers[i] == 1) {
        scores1 <- c(NA, scores1)
        scores2 <- c(NA, scores2)
      } else {
        if(outliers[i] < length(scores1)) {
          scores1 <- c(scores1[1:(outliers[i]-1)], NA, scores1[outliers[i]:length(scores1)])
          scores2 <- c(scores2[1:(outliers[i]-1)], NA, scores2[outliers[i]:length(scores2)])
        } else {
          scores1 <- c(scores1, NA)
          scores2 <- c(scores2, NA)
        }
      }
    }
  } else {
    cat('Not enough genes to remove even 1 loci, use different outlier_fraction', '\n')
    scores1 <- PCA$scores[,1]
    scores2 <- PCA$scores[,2]
  }
} else {
  scores1 <- PCA$scores[,1]
  scores2 <- PCA$scores[,2]
}

variables$PC_1 <- scores1
variables$PC_2 <- scores2
if(any(is.na(variables))) variables <- variables[which(complete.cases(variables)),]

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
  if(as.numeric(unlist(strsplit(PC_rate, '_'))[2]) == 1) {
    PC_usefulness <- 2
  } else {
    PC_usefulness <- 1
  }
  
  loadings_usefulness = loadings(PCA)[][,PC_usefulness]
  if(loadings_usefulness[length(loadings_usefulness)] < 0 && loadings_usefulness[length(loadings_usefulness)-1] < 0) {
    if(any(loadings_usefulness[1:(length(loadings_usefulness)-3)] > 0)) {
      direction <- 'clear'
      descending <- F
    } else {
      direction <- 'unclear'
    }
  } else {
    if(loadings_usefulness[length(loadings_usefulness)] > 0 && loadings_usefulness[length(loadings_usefulness)-1] > 0) {
      if(any(loadings_usefulness[1:(length(loadings_usefulness)-3)] < 0)) {
        direction <- 'clear'
        descending <- T
      }
      direction <- 'unclear'
    } else {
      direction <- 'unclear'
    }
  }
} else {
  direction <- 'unclear'
}

#D) Sort & Subsample------------------------------------------------------------------------------
if(direction == 'unclear') {
  cat('It is unclear how to sort the data.', '\n', 'You can the check loadings and decide manually how to proceed.', '\n',  
      'It is also possible that this approach might not be suitable for your dataset.', '\n', 
      'If this is the case, choose a different gene property with which to sort.', '\n')
}

if(grepl('maybe', PC_rate)) {
  cat('There seems to be some ambiguity as to the identity of the axes.', '\n', 
  'Proceed with caution and check PCA loadings to see if sorting is appropriate')
}

if(direction == 'clear') {
  usefulness_col <- grep(PC_usefulness, colnames(variables))
  
  #sort by usefulness
  variables_sorted <- variables[order(variables[,usefulness_col], decreasing = descending),]
}

###WARNING: uncomment and omdify the following lines if you would like to sort
###and subsample by a different property, for example occupancy or RF similarity
#variables_sorted = variables[order(variables[,'occupancy'], decreasing = T),]
#variables_sorted = variables[order(variables[,'robinson_sim'], decreasing = T),]

#sort entire dataset according to the sorting order imposed
positions = c()
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

#subsample (or not if n_genes was left as = 0)
sorted_names <- sorted_names[1:n_genes]
sorted_partitions <- sorted_partitions[1:n_genes,]
sorted_data <- sorted_data[,1:sorted_partitions$End[n_genes]]
sorted_trees <- sorted_trees[1:n_genes]

write.phyDat(as.phyDat(sorted_data), file = paste0(getwd(), '/sorted_alignment_', n_genes, 'genes.fa'), format = 'fasta')
partitions_tosave = paste0(sorted_names, variables_sorted$genes[1:n_genes], ' = ', sorted_partitions$Start, '-', sorted_partitions$End)
write(partitions_tosave, file = paste0(getwd(), '/sorted_alignment_', n_genes, 'genes.txt'))
write.tree(sorted_trees, file = paste0(getwd(), '/sorted_trees_', n_genes, 'genes.txt'))

#E) Optional: visualize some sorting results
variables_to_plot <- data.frame(gene = rep(variables_sorted$genes, ncol(PCA$loadings)),
                               value = c(as.matrix(variables_sorted[,variables_to_use])),
                               property = rep(colnames(variables_sorted[,variables_to_use]), each = nrow(variables_sorted)),
                               pos = rep(1:nrow(variables_sorted), ncol(PCA$loadings)))

order_properties <- as.character(unique(variables_to_plot$property))

if(nrow(partitions) == nrow(sorted_partitions)) {
  ggplot(variables_to_plot, aes(x = pos, y = value, color = factor(property, levels = order_properties))) + 
    geom_smooth() + theme_bw() + theme(legend.position = "none") + 
    facet_wrap(vars(factor(property, levels = order_properties)), scales = 'free', nrow = 2)
} else {
  ggplot(variables_to_plot, aes(x = pos, y = value, color = factor(property, levels = order_properties))) + 
    geom_smooth() + theme_bw() + theme(legend.position = "none") + 
    facet_wrap(vars(factor(property, levels = order_properties)), scales = 'free', nrow = 2) +
    geom_vline(xintercept = n_genes, linetype = 'dashed')
}

#Other funtions-----------------------------------------------------------------------------------
`%not in%` <- function(x, table) is.na(match(x, table, nomatch=NA_integer_))
#function to count invariant sites
inv <- function(x) {
  pattern <- unique(x)
  if(any(pattern %in% c('?'))) pattern <- pattern[-which(pattern == '?')]
  if(any(pattern %in% c('-'))) pattern <- pattern[-which(pattern == '-')]
  
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
    missing = which(names(unlist(x)) == '-')
    x <- x[-missing]
  }
  if('?' %in% names(unlist(x))) {
    missing <- which(names(unlist(x)) == '?')
    x <- x[-missing]
  }
  return(x)
}