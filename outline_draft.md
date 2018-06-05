Bacterial pan-genome analysis with population graphs

Intro

1. The biology - what do we expect bacterial genomes to look like?
   - In bacteria, variation occurs in 3 ways: mutation, homologous recombination non-homologous recombination
   - For some species (e.g. E coli) the mutations occur as often as  homolgous recombinations
   - Mutations result in a closely related sequence
   - Homolgous recombinations give rise to gene sequences which are mosaics of previously seen gene sequences
   - Non-homologous recombination results in a genome with a large degree of structural and content variability
2. This underlying biology and the availability of sequenced genomes context ideal for pangenome graphs which explicitly infer new sequences as a mosaic of a panel of reference sequences
3. What do we want to use genome sequence data for?
   - Understanding genetic evolution of a species
   - Investigating why different strains of a species may have different phenotypes
   - In a clinical setting, for outbreak tracking and AMR prediction
4. What to current methods allow us to infer from sequence data and what do are their limitations?
   - Gene presence/absence (AMR prediction, investigating recombination) e.g. 
     For Illumina data, downstream analysis can then allow variation detection in specific genes of interest (phenotype deconvolution, AMR classification)
   - Variation detection in the core genome from Illumina (outbreak tracking)
   - Whole genome MSA (doesn't scale, not robust to rearrangements)
   - wgMLST - does allow classification as a series of barcodes, but doesn't allow relatedness between barcodes to be easily defined
   - clairvoyant only tool for snv typing from nanopore, but designed for human and requires close reference
   - None provide a framework to combine these and to describe variation in both the core and accessory genome
   - None allow variation detection from nanopore in bacteria
5. What do we do?
   - We introduce a novel graph structure and inference method which allows simul-
     taneous  identification  of  pangenome  elements  which  are  present  and  genotyping
     within them
   - This enables SNP and indel calling in all pangenome genes
   - We handle long read sequence data, and  is  backward  compatible  with  illumina.  
   - We  are  able  to make  accurate  calls  from  nanopore  data  (cardio) and note that this is not species specific (tb)
   - We  add  resolution  in  the  context  of  an  E.  coli  outbreak  by  considering  the  whole
     pangenome.
   - We demonstrate usefullness in a clinical setting by allowing comparison of isolates which are environmentally linked but may not be genetically related (nicole illumina prospective).

Results

1. Algorithms
   - Novel pan genome reference structure
   - Indexing PRGs with minimizing kmers
   - Inferring sequences as a mosaic through graph
   - Describing variation by comparison with an optimal reference
   - Novel variation detection
2. Datasets
   - Ecoli K12 - proof of principle, sample in graph, equivalent to simulations
     - Prove sequences do lie in graph
     - Given that sequences lie in graph, is like running simulations
     - Note don't need novel caller for this section
     - Proportion of genes correctly identified as present, number of false positives by coverage
     - For SNP sites in true positive gene vcfs, what proportion genotyped correctly by coverage
     - For indel sites in true positive gene vcfs, what proportion genotyped correctly by coverage
     - For all sites in true positive gene vcfs, what proportion genotyped correctly by coverage
     - For inferred gene sequences, how many bases away from reference genome is it?
     - Overall inferred gene sequence accuracy
     - Method:
       1. Only work with annotated genes which we have checked exist exactly in the graph
       2. For each 10x coverage, run pandora
       3. Plot number of annotated genes found and number of genes which did not map back to reference vs coverage
       4. Plot % GT correct for snps, indels, all vs coverage
       5. Supplementary plot for frequency distribution of number of mismatch bases in genes
   - TB - accuracy from nanopore/novel detection
   - Cardio - comparison with other tools
     - To describe outbreak, want a multisample vcf showing how samples differ and a tree
     - Nanopore comparison method: Starting with the same read fasta, assemble with canu, polish with nanopolish and use minimap2/bwa to output a vcf with respect to (a) an ecoli reference, (b) the pacbio sequenced outbreak reference. Combine these accross samples and build a tree from them
     - For pandora, force vcf output to have reference from (a), (b) and try with best reference
     - Illumina comparison method: Cortex the samples to get a multivcf and turn into a tree

Discussion

- Implications
- Fundamental limitations
- Non-fundamental limitations

Conclusions

Methods



Supplemental Material

1. Performance



Thoughts about comparitors

1. PanSeq takes genomes or draft genomes and identifies regions which are core or accessory to them. Within the core it finds SNPs, and for all samples outputs presence/absence of all regions
   +ve identifies regions accross pangenome
   -ve only identifies SNPs
   -ve only variant calls in core
   -ve takes draft genomes (but this means could handle nanopore?)
2. Harvest (ParSNP) does core genome alignment and finds all variants in this core genome, with visuals from gingr
   +ve identifies regions accross pangenome
   +ve identifies SNPs indels and structural variants
   -ve only variation in core
   -ve takes draft genomes
   -ve requires genomes to be very closely related e.g. >97% id
3. Cortex allows variant calls between different samples (e.g. outbryk) everywhere
   +ve no assembly
   +ve all read data
   -ve complexity
   -ve no nanopore error handling
4. JR assembled read data (spades, hybridspades, canu) then mapped using their pipeline? The MSA output of this passed to IQTree and ClonalFrameML for phylogenies.
5. ARIBA?
6. GROOT - variation aware reference graphs for resistome typing
   -ve doesn't look like it does anything other than type

Comments from Zam which need to be integrated by Rachel:

- There is a market for understanding bacterial genomes in general, not just in a healthcare setting. e.g. understanding why phenotypes are different, or how species evolve.
- It is self evident that looking at SNPs accross ALL genes is better and in this context we make a gain as a result
- Main demonstrated use case - the prospective illumina dataset of Nicole's would be ideal provided that Nicole or someone else who knows about the dataset could look over the results with us. Want to see if clustering based on pandora SNPs is closed to epi data than alternative method
- Demonstration of accuracy - cardio may be better than TB as have more nanopore samples with truth
- Definitely need to include TB accuracy if it turns out that there is a different error bias
- MSA doesn't scall and not robust to rearrangements, difficult to do
- wgMLST more like what we do, but it is a pain to relate all the different barcodes
- We can in theory with multisample output info about there being 2 refs with 80% samples close to the first
- We should aim to output a fastq with low qual where we are unsure, rather than a fasta of gene sequences
- When using cardio to demonstrate can make accurate cakks, compare with samtoools or snappy, not clockwork
