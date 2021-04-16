# genesortR
Sorting and subsampling of phylogenomic datasets using a multivariate method to quantify phylogenetic usefulness.

## Description
This R script estimates seven gene properties commonly used to characterize the information content of loci in phylogenomic datasets: four sources of systematic bias (average pairwise patristic distance, compositional heterogeneity, level of saturation, and root-to-tip variance), two proxies for phylogenetic signal (Robinson-Foulds similarity to a target topology and average bootstrap support), as well as the proportion of variable sites. Some of these properties are estimated directly on sequence data, others on their corresponding topologies. Instead of directly optimizing these properties, a principal component analysis (PCA) is used to find an axis of phylogenetic usefulness along which proxies for signal increase while sources of bias decrease.

This approach was first used for phylogenomic subsampling by Mongiardino Koch & Thompson (2020), and then found to be able to discover an axis of phylogenetic usefulness across a sample of 18 diverse phylogenomic datasets by Mongiardino Koch (2021). The method was also found to recover more accurate topologies than most other commonly-used methods for phylogenomic subsampling, including those based on rates of evolution (targeting either low, intermediate or high rates), as well as methods that seek to minimize potential sources of systematic bias (e.g., targeting the most clock-like, least saturated, and least compositionally heterogeneous loci). It therefore provides a well-founded alternative to reduce the size of phylogenomic datasets, allowing the use of more complex inference methods (full coalescent methods, site-heterogeneous models, total-evidence dating, etc.), as well as allowing phylogenetic hypotheses to be tested with smaller and better curated datasets.

## Parameters
The script requires some editing to define a number of input parameters:
* Set the woring directory to a desired location, ideally where input files are located. This is also where output files will be written.
* Provide names for the four data files that need to be loaded:
  1. Data matrix: needs to be in FASTA format. Only '?', '-' and 'X' are taken to represent missing data. If missing data is represented with other symbols, consider editing the data before running the script. 
  2. Partition file: in format 'gene_name = start-end'.
  3. Gene trees: a single tree file in newick format with all gene trees (with node support values).
  4. Species tree: a target topology used to calculate Robinson-Foulds similarity (also in newick format). If some relationships between the sampled taxa are contentious, these can be collapsed so as not to favor any specific resolution. IMPORTANT: This species tree has to be rooted using outgroups.
* ```type```: the data type, either 'AA' or 'DNA'.
* ```ingroup``` (OPTIONAL): The names of two terminals that bracket the ingroup clade. This is recommended and serves to purposes: it allows gene trees to be rooted (which is important for the calculation of root-to-tip distances for example), and it allows outgroups to be discarded before measuring gene properties. If this is not desired, then this parameter can be left empty (i.e., unmodified) and all terminals will be used in the estimation of properties. If terminal names are provided, gene trees will be rooted at the node of the most-recent common ancestor (MRCA) of the ingroup.
* ```threshold```: A threshold of number of ingroup taxa (i.e., a level of occupancy within the ingroup) to even consider the loci. If left as default (i.e., 'auto'), then this is set to 10% of the ingroup, otherwise it can be modified to a specific number of taxa. Loci with less than ```threshold``` will be removed from the data. If this is not desired set to 0.
* ```remove_outliers```: A logical value stating whether to activate the removal of outlier loci (recommended if the dataset is composed of hundreds of loci or more). This will discard a fraction of loci (defined using ```outlier_fraction```) that exhibit the largest Mahalanobis distance in PC space, and which may suffer from problems in orthology inference, alignment, etc. Even if no issue can be detected with these loci, their pattern of correlation between gene properties deviates markedly from that of all others, indicating unusual evolutionary histories (e.g. driven by strong selective pressures) and it is unlikely that small datasets will benefit from their inclusion. The PCA is then repeated after their exclusion.
* ```n_genes```: A final number of genes to retain. If left as defualt (i.e., 'all') then the data is sorted but not subsampled, and all loci are saved to output files. If the objective is to subsample the dataset, modify to a desired number of loci in the final dataset.
* ```topological_similarity```: A logical value determining whether to incorporate Robinson-Foulds (RF) similarity in the PCA or not. This is highly recommended, as in the absence of RF similarity the only proxy for phylogenetic signal left is the average bootstrap support. I have found that, across several datasets, this can lead to favoring the selection of loci that evolve faster than necessary (and contain a higher prevalence of some issues such as compositional and rate heterogeneity). The addition of RF similarity helps balance this effect. However, this requires using a species tree that is considered "true", i.e., it is used as a target against which the topologies of gene trees are evaluated. In the presence of extensive phylogenetic uncertainty, this can be considered undesirable, and potentially even mislead the selection of loci. While even in the case of a highly uncertain phylogeny I would still recommend the inclusion of RF similarity, possibly by providing a relatively unresolved topology with uncertain nodes collapsed, ```topological_similarity``` can be set to ```FALSE``` and RF will not be used to constrain the phylogenetic usefulness axis.<br/><br/>
**WARNING: Taxon names need to match across all files, and loci need to be ordered in the same way in the alignment, partition and gene tree files.**

## Estimated gene properties
Seven gene properties are calculated in order to derive PC axes, inlcuding four sources of systematic bias (1-4 below) and two proxies for phylogenetic signal (5-6 below). In every case, sources of bias and proxies for signal are defined in such a way that their minimization and mazimization (respectively) is desireable.

1. Root-to-tip variance: This metric provides an estimate of the clocklikeness of the evolutionary process for a loci, andh has been employed routinely to select genes for divergence time estimation (e.g., Smith et al. 2018). If an ```ingroup``` is defined by providing terminal names that bracket the clade, then gene trees will be rooted at the node of their MRCA, and distances will be calculated from this root all terminals within the ingroup. Distances to outgroup taxa will not be used, even when these form part of the ingroup clade in a given gene tree. If an ingroupis not defined, the tree will be rooted at its mid-point and distances will be estimated from this root to all terminals.
2. Average patristic distances: This value is estimated as the average of all pair-wise patristic distances (i.e., the sum of the lengths of branches that link two nodes in the gene tree) between terminals. Higher values are considered to be conducive to long-branch attraction artifacts (Struck 2014, Kocot et al. 2017).
3. Level of saturation: Estimated as the complement of the regression slope of pair-wise patristic distances (see above) against p-distances (the proportion of sites in the alignment for which two sequences differ) (Nosenko et al. 2013). Although negative correlations can occur, these are taken to happen due to some level of random error in a completely saturated loci, and the value is replaced with a 1.
4. Compositional heterogeneity: Variation in the use of amino acids among different branches of a tree is a generally unaccounted source of non-phylogenetic signal, and a major source of error for some clades (Nesnidal et al. 2013). This heterogeneity is estimated here using the relative composition frequency variability (RCFV), as defined by Nesnidal et al. (2010). Note that this variable is not included in the PCA if data type is set to 'DNA'.
5. Robinson-Foulds (RF) similarity: This variable (which corresponds to the complement of the RF distance, Robinson & Foulds 1981) is used as an estimate of gene tree error. Although this can suffer from issues of circularity, it should be noted that the species tree provided (which is used as the target topology for the calculation of RF similarity scores) can contain any number of unresolved nodes depicting uncertainty in the resolution of clades.
6. Average bootstrap support: The average degree of support in the molecular data for the nodes in the gene tree.
7. Proportion of variable sites: This variable has been variously consider an estimate of information content, a proxy for evolutionary rate, or an amount of phylogenetic signal in the alignment. Its degree of correlation with other gene properties across multiple datasets supports the first two interpretations.

Besides these seven properties directly employed for the discovery of a phylogenetic usefulness axis, six other properties are also estimated although not directly employed. These include: A. The amount of missing data; B. The level of occupancy; C. The length of the alignment; D. The total tree length (i.e. the sum of all of its branches); E. An estimate of the rate of evolution (tree length divided by the number of taxa; Telford et al. 2014); and F) The treeness of a gene tree (i.e., the proportion of branch lengths in internal branches).
Any of these thirteen properties can be directly employed for the sorting of datasets with a minimial degree of editing to the R script. Instructions on how to do this are provided as comments in section 'D) Sort & Subsample' of the script.

## Output
The main outputs that are directly saved to file are the sorted alignment, partition file and gene tree file. If a desired number of genes is specified using the ```n_genes``` parameter, then these datasets will be subsampled as well, allowing for both concatenation and coalescent-based phylogenetic inference to be repeated with a smaller dataset. The estimated properties for both the outlier loci (if outlier filtering is activated), and the full sorted dataset are also output as csv files (properties_outliers.csv and properties_sorted_dataset.csv, respectively). This can help explore the reasons why loci were removed (and potentially fix some of their issues if these are evident), as well as check the correlation among different gene properties and their distribution along the ordering imposed by the script. Additionally, a plot (Fig. 1) is generated showing how the underlying gene properties vary according to the order in which genes are placed based on their phylogenetic usefulness. This plot is also saved to file.<br/><br/><br/>
![sorting_example](https://github.com/mongiardino/genesortR/blob/main/images/sorting_example.jpeg)
**Fig. 1:** Value of the seven gene properties against the order in which loci are sorted according to their phylogenetic usefulness. Loci are from a sea urchin transcriptomic dataset (2,356 loci) of Mongiardino Koch & Thompson (2020). Regression lines correspond to generalized additive models (GAM).

## Author
Nicolás Mongiardino Koch. Department of Earth & Planetary Sciences, Yale University.

Citation: Mongiardino Koch N. 2021. Phylogenomic subsampling and the search for phylogenetically reliable loci. bioRxiv 2021.02.13.431075. https://doi.org/10.1101/2021.02.13.431075.

I greatly appreciate the efforts of Mansa Srivastav and Kevin Kocot in helping debug this script.

## References
Kocot KM, Struck TH, Merkel J, Waits DS, Todt C, Brannock PM, Weese DA, Cannon JT, Moroz LL, Lieb B, Halanych KM. 2017. Phylogenomics of Lophotrochozoa with consideration of systematic error. Systematic Biology 66:256-282.

Mongiardino Koch N, Thompson JR. 2020. A total-evidence dated phylogeny of Echinoidea combining phylogenomic and paleontological data. Systematic Biology syaa069, https://doi.org/10.1093/sysbio/syaa069.

Mongiardino Koch N. 2021. Phylogenomic subsampling and the search for phylogenetically reliable loci. bioRxiv 2021.02.13.431075. https://doi.org/10.1101/2021.02.13.431075.

Nesnidal MP, Helmkampf M, Bruchhaus I, Hausdorf B. 2010. Compositional heterogeneity and phylogenomic inference of metazoan relationships. Molecular Biology and Evolution 27:2095-2104.

Nesnidal MP, Helmkampf M, Meyer A, Witek A, Bruchhaus I, Ebersberger I, Hankeln T, Lieb B, Struck TH, Hausdorf B. 2013. New phylogenomic data support the monophyly of Lophophorata and an Ectoproct–Phoronid clade and indicate that Polyzoa and Kryptrochozoa are caused by systematic bias. BMC Evolutionary Biology 13:253.

Nosenko T, Schreiber F, Adamska M, Adamski M, Eitel M, Hammel J, Maldonado M, Müller WE, Nickel M, Schierwater B, Vacelet J. 2013. Deep metazoan phylogeny: when different genes tell different stories. Molecular Phylogenetics and Evolution 67:223-233.

Robinson DF, Foulds LR. 1981. Comparison of phylogenetic trees. Mathematical Biosciences 53:131-730.

Smith SA, Brown JW, Walker JF. 2018. So many genes, so little time: A practical approach to divergence-time estimation in the genomic era. PloS One 13:e0197433.

Struck TH. 2014. TreSpEx – Detection of misleading signal in phylogenetic reconstructions based on tree information. Evolutionary Bioinformatics 10:51-67.

Telford MJ, Lowe CJ, Cameron CB, Ortega-Martinez O, Aronowicz J, Oliveri P, Copley RR. 2014. Phylogenomic analysis of echinoderm class relationships supports Asterozoa. Proceedings of the Royal Society B 281:20140479.
