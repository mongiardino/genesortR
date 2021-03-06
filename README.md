# genesortR
Sorting and subsampling of phylogenomic datasets using a multivariate method to quantify phylogenetic usefulness.

## Description
This R script estimates seven gene properties commonly used to characterize the information content of loci in phylogenomic datasets: four sources of systematic bias (average pairwise patristic distance, compositional heterogeneity, level of saturation, and root-to-tip variance), two proxies for phylogenetic signal (Robinson-Foulds similarity to a target topology and average bootstrap support), as well as the proportion of variable sites. Some of these properties are estimated directly on sequence data, others on their corresponding topologies. Instead of directly employing these to reduce the size of datasets, a principal component analysis is used to find an axis of phylogenetic usefulness along which proxies for signal increase while sources of bias decrease. This approach was first used for phylogenomic subsampling by Mongiardino Koch & Thompson (2020), and then found to be an axis prevalent across a sample for 18 diverse phylogenomic datasets by Mongiardino Koch (2021).

![model_comparison](https://github.com/mongiardino/genesortR/tree/main/images/sorted_example.jpg)
**Fig. 1:** Value of the seven gene properties against the order of loci sorted by phylogenetic usefulness.
