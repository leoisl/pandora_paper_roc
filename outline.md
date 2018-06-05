__Nucleotide-resolution bacterial pan-genomics with reference  graphs__


__Intro__

Want to access all variation, exists in layers

__Results__

1. Construction of PRG, (w,k)-minimizers
2. False positive/negative rates for gene presence detection as function of gene length and error rate
   (I assume the right thing to do is ignore genes below some length?)
3. Error rate and profile for inferred mosaic
4. Adding nanopolish
5. report variation with per-gene best reference


__Application__

1. Cardio-dataset: how well does it do on sequence inference and gene inference. What does tree look like?
2. Ten M. tuberculosis genomes - can we get a good set of SNPs and indels? What depth of coverage do we need?

__Superiority over existing methods__

1. Avoids using MSA, avoids using core.
2. Compare whole genome assembly and MSA/nanopolish.
3. Remove nanopore systematic errors by using prior (graph)


__Things I propose ignoring__

1. gene order inference
2. Application: AMR gene on plasmid context? V hard - see here http://aac.asm.org/content/60/6/3767/F2.large.jpg
