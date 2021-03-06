# genesortR
Sorting and subsampling of phylogenomic datasets using a multivariate method to quantify phylogenetic usefulness.

## Description
This R script estimates seven gene properties commonly used to characterize the information content of loci in phylogenomic datasets: four sources of systematic bias (average pairwise patristic distance, compositional heterogeneity, level of saturation, and root-to-tip variance), two proxies for phylogenetic signal (Robinson-Foulds similarity to a target topology and average bootstrap support), as well as the proportion of variable sites. Some of these properties are estimated directly on sequence data, others on their corresponding topologies. Instead of directly optimizing these properties, a principal component analysis (PCA) is used to find an axis of phylogenetic usefulness along which proxies for signal increase while sources of bias decrease. This approach was first used for phylogenomic subsampling by Mongiardino Koch & Thompson (2020), and then found to be an axis prevalent across a sample of 18 diverse phylogenomic datasets by Mongiardino Koch (2021). This approach was found to recover more accurate topologies than most other commonly-used methods for phylogenomic subsampling, include those based on rates of evolution (targeting either low, intermediate or high rates), as well as methods that seek to minimize potential sources of systematic bias (e.g., targeting the most clock-like, least saturated, and least compositionally heterogeneous loci). The method thus provides a well-founded alternative to reduce the size of phylogenomic datasets, allowing the use of more complex inference methods (full coalescent methods, site-heterogeneous models, total-evidence dating, etc.), as well as allowing phylogenetic hypotheses to be tested with smaller, better curated datasets.

## Parameters
The script requires some editing to define a number of input parameters:
* Set the woring directory to a desired location, ideally where input files are located. This is also where output files will be written.
* Provide names for the four data files that need to be loaded:
  1. Data matrix: needs to be in FASTA format.
  2. Partition file: in format 'gene_name = start-end'.
  3. Gene trees: a single tree file in newick format with all gene trees (with node support values).
  4. Species tree: a target topology used to calculate Robinson-Foulds similarity (also in newick format). If some relationships between the sampled taxa are contentious, these can be collapsed so as not to favor any specific resolution.
* ```type```: the data type, either 'AA' or 'DNA'.
* ```ingroup``` (OPTIONAL): The names of two terminals that bracket the ingroup clade. This is recommended and serves to purposes: it allows gene trees to be rooted (which is important for the calculation of root-to-tip distances for example), and it allows outgroups to be discarded before measuring gene properties. If this is not desired, then this parameter can be left empty (i.e., unmodified) and all terminals will be used in the estimation of properties.
* ```threshold```: A threshold of number of ingroup taxa (i.e., a level of occupancy within the ingroup) to even consider the loci. If left as default (i.e., ```auto```), then this is set to 10% of the ingroup, otherwise it can be modified to a specific number of taxa.
* ```remove_outliers```: A logical value stating whether to activate the removal of outlier loci (recommended if the dataset is composed of hundreds of loci or more). This will discard a fraction of loci (defined using ```outlier_fraction```) that exhibit the largest Mahalanobis distance in PC space, and which may suffer from problems in orthology inference, alignment, etc. Even if no issue can be detected with these loci, their pattern of correlation between gene properties deviates markedly from that of all others, indicating unusual evolutionary histories (e.g. driven by strong selective pressures) and it is unlikely that small datasets will benefit from their inclusion. The PCA is then repeated after their exclusion.
* ```n_genes```: A final number of genes to retain. If left as defualt (i.e., ```all```) then the data is sorted but not subsampled, and all loci are saved to output files. If the objective is to subsample the dataset, modify to a desired number of loci in the final dataset.<br/>
**Note that taxon names need to match across all files, and that loci need to be ordered in the same way in the alignment and gene tree files.**
<br/><br/><br/>

![sorting_example](https://github.com/mongiardino/genesortR/blob/main/images/sorting_example.jpeg)
**Fig. 1:** Value of the seven gene properties against the order of loci sorted by phylogenetic usefulness.
