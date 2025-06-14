# ANCOM-implementation
An implementation in R of ANCOM(Mandal et al. 2015), using a custom method for finding a Reference taxon for computation

Please not: this implementation is designed for populations(K)=2 and subjects per population(j)=1.

Upon taxonomic characterisation of 16s rRNA sequencing samples, we obtain Operational Taxonomic Units(OTUs) or Amplicon Sequence Variants(ASVs) which represent a unit of a microbial taxa.

These units usually represent Relative abundance of taxons in a community. While Absolute abundance is count data which can be applied directly to parametric statistical methods, Relative abundance are counts out of a total whole. This is because Amplicons produced DNA sequencing are automatically normalised by the sequencing depth or number of total reads produced.
This kind of data is capped upto a maximum value and is present in a simplex space. If ASV of one species of organism increases, it corresponds to a proportionate decrease in ASVs of one or multiple other species. Another example of this data type is percentage data(out of 100).

Relative abundance cannot be normally distributed, and is thus not directly amenable to statistical tests such as linear regression, T-tests, ANOVA and others.
Also, the total number of organisms in these samples may not be known, along with the sampling fractions.
This data usually needs to transformed in some manner before being subjected to statistical tests. 

In an experiment where the abuhndace of a specific ASV in sample 1 needs to be compared with the abundance of the same ASV in sample 2, a direct T-test cannot be used directly.
This project uses the methods of Mandal et al.(https://www.tandfonline.com/doi/full/10.3402/mehd.v26.27663) to enable the discovery of taxons which are significantly different between two samples by using a log transform method on Relative abundance of ASVs.

Since the reference NACOM paper does not seem to describe a method for findinga reference taxon that does not vary in Absolute abundance between the two samples,
a custom mathematical formula was developed and utilised here, which is as follows:

(Log(E(γr1))-Log(E(γr2)))^2

Where:

γr1 represents the relative abundance of taxon 1 in sample 1, 
γr2 represents the relative abundance of taxon 1 in sample 2, 
E(x) represents the average or expected value of x.

Here log ratio of relative abundances were taken to allow for Absolute abundance based comparision, which was then squared to account for bidirectional differences. The taxon outputting the lowest value from this formula was taken as the reference taxon.

# Inputs

Inputs need to be two tables, each from a specific time point, containing the relative abundances for each observed taxon accross all samples. The tables need to specify the taxon names in the columns and the sample names in the rows, as illustrated:

samples  Taxon1 Taxon2 Taxon3
Sample1  0      2      4
Sample2  7      2      9

Please not that the union of taxons from all samples needs to be present in the table. Missing values can be filled with zeoes. Also the order of the taxon names and the samples needs to be the same between both tables.

For detailed methods please refer to the ANCOM paper: https://www.tandfonline.com/doi/full/10.3402/mehd.v26.27663.

Disclaimer: This method uses the following assumptions, as far as the author was able to glean:
1. PCR amplication of rRNA is not subjected to sequencing bias.
2. All microbes are uniformly distributed accross the samples.
3. Atleast one taxon does not vary in terms of Absolute counts between the two samples(populations).
4. Difference in cell counts between the samples is not the same for every taxon.
