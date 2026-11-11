# Supplementary methods: restricted mitochondrial tRNA-Phe(GAA) locus

## Scope and selection

A targeted follow-up identified 16 vertebrate mitochondrial tRNA-Phe(GAA) candidates within the RF00005 seed and checked each candidate against an accession-versioned complete mitochondrial genome record. The accepted set is restricted to this one mitochondrial locus; it is not a claim that all RF00005 sequences, or all tRNA-Phe genes, are orthologous.

## Validation criteria

Each candidate was retained only when the available record-level evidence supported the same locus assignment:

1. The ungapped Rfam Stockholm sequence matched the corresponding GenBank tRNA-Phe feature sequence.
2. The mapped Rfam RF model positions 33–35, corresponding to Stockholm columns 44–46, contained the GAA anticodon.
3. The complete mitochondrial genome contained 22 annotated tRNA features and exactly one tRNA-Phe feature.
4. The local gene order was tRNA-Phe → 12S rRNA → tRNA-Val → 16S rRNA.

All 16 candidates met the sequence-feature, complete-genome, copy-number, and gene-order checks. Three records carried explicit GenBank anticodon qualifiers; all three agreed with the alignment-coordinate mapping. The [acceptance record](provenance/trna_phe_orthology_acceptance.json) summarizes these counts and the restricted decision. The [validation table](tables/trna_phe_locus_validation.csv) lists each organism and its Rfam and complete-mitogenome accession.version identifiers, sequence identity and coverage, tRNA feature count, tRNA-Phe copy count, and gene-order result.

## Descriptive analysis and sensitivity

The subset contains 42 paired and 29 explicitly unpaired alignment columns. Across 21 consensus pairs, 324 of 336 possible sequence-by-pair observations were valid; 305 of those 324 observations were AU, UA, GC, CG, GU, or UG. The paired-minus-unpaired conservation contrast is sensitive to the valid-column threshold: it changes from +0.036 at the full set to −0.006 when the minimum valid fraction is raised to 0.90. The one-Bovidae-representative analysis also changes the contrast across the nine possible choices. This is a taxon-composition sensitivity analysis, not a phylogenetic correction.

## Data sources and limitations

The accession.version identifiers in the validation table link to public NCBI nucleotide records. The underlying NCBI sequence records are not copied into this repository. Feature and anticodon-location interpretation follows the [DDBJ feature-table specification](https://www.ddbj.nig.ac.jp/ddbj/feature-table-e.html) and [ENA anticodon qualifier reference](https://www.ebi.ac.uk/ena/WebFeat/qualifiers/anticodon.html).

No species-tree-aware evolutionary analysis is included. Conservation and base-pair composition alone do not establish compensatory evolution, causation, or functional constraint. The subset result must not be used to replace the broad RF00005 family summary.
