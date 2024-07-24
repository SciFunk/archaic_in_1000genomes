# Archaic allele counts in 1000 Genomes Populations

The script sprime_conversion.py takes the archaic allele predictions from [Sprime](https://doi.org/10.1016/j.cell.2018.02.031) and converts them to a bed file for use with vcftools. vcftools can then be used to extract the archaic sites from [1000 Genomes Project data](https://www.internationalgenome.org/data-portal/data-collection/phase-3).

The script sprime_allele_counts.py calculates the archaic allele frequency for each of the 26 populations in the 1000 Genomes Project.
