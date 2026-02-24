"""
Variant scoring functions using AlphaGenome.
"""

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers
from . import _get_dna_model


def score_variant(chrom, pos, ref, alt):
    """
    Score a single variant using AlphaGenome.
    
    Args:
        chrom: Chromosome (e.g., 'chr17')
        pos: Position (1-based, hg38)
        ref: Reference allele
        alt: Alternate allele
    
    Returns:
        DataFrame with columns:
            - gene_name: affected gene
            - output_type: RNA_SEQ, ATAC, etc
            - ontology_curie: tissue identifier
            - biosample_name: tissue name
            - raw_score: effect size (log fold change for RNA_SEQ)
            - quantile_score: percentile ranking
    """
    dna_model = _get_dna_model()
    
    # Define the variant
    variant = genome.Variant(
        chromosome=chrom,
        position=pos,
        reference_bases=ref,
        alternate_bases=alt,
    )

    # Build 1MB window around the variant
    interval = variant.reference_interval.resize(
        dna_client.SUPPORTED_SEQUENCE_LENGTHS['SEQUENCE_LENGTH_1MB']
    )

    # Score using all recommended scorers (RNA_SEQ, ATAC, CAGE, etc)
    variant_scores = dna_model.score_variant(
        interval=interval,
        variant=variant,
        variant_scorers=list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values()),
    )

    # Convert to clean dataframe
    df_scores = variant_scorers.tidy_scores(variant_scores)

    return df_scores
