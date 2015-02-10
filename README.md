# ParTI

Biological research increasingly depends on interpreting large datasets in high dimensional space. The main approaches for analyzing such data are dimensionality reduction techniques such as principal component analysis, and clustering techniques that split data points into groups. Recently, a theoretical advance (Shoval et al., Science 2012) has suggested a complementary way to understand large biological datasets. This approach is based on Pareto optimality of organisms with respect to multiple evolutionary tasks. It predicts that cells or organisms that need to perform multiple tasks have phenotypes that fall on low dimensional polytopes such as lines, triangles, tetrahedrons; phenotypes optimal for each task – called archetypes – are at the vertices of these polytopes.

The present software package implements the Pareto Task Inference (ParTI) method to analyze biological data in light of this theory. ParTI finds the best fit polytopes and the corresponding archetypes, and indicates which features of the data are enriched near each archetype. The present approach is less sensitive to the density of sampling of the data space than clustering approaches, because it relies on the overall shape of the data and its pointy vertices; enrichment of relevant biological features is higher near the archetypes than near the centers of clusters obtained by standard methods.

The package comes with two example datasets and template scripts to analyse them: a dataset of 2000 tumor breast cancer gene expression dataset by Curtis et al., Nature 486, 345 (2012), and 63 mouse tissues by Lattin et al., Immunome Res 4, 5 (2008).

The article describing the approach can be found here:
http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.3254.html

The documentation of the package can be found here:
http://www.weizmann.ac.il/mcb/UriAlon/download/ParTI
