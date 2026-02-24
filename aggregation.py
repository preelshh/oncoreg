"""
Functions for aggregating variant effects into gene-level and patient-level metrics.
"""

import pandas as pd
from .data import get_tissue_ontology


def compute_gene_impact(list_of_df_scores, cancer_type='breast'):
    """
    Aggregate variant effects across multiple variants to compute per-gene impact.
    
    Args:
        list_of_df_scores: List of DataFrames from score_variant()
        cancer_type: Cancer type ('breast', 'lung', 'colon', etc)
    
    Returns:
        DataFrame with columns:
            - gene_name: gene symbol
            - gene_impact: mean absolute effect across all variants affecting this gene
    """
    # Get tissue ontology for this cancer type
    ontology = get_tissue_ontology(cancer_type)
    
    # Concatenate all variant scores
    all_scores = pd.concat(list_of_df_scores, ignore_index=True)
    
    # Filter to target tissue and RNA_SEQ only
    filtered = all_scores[
        (all_scores['ontology_curie'] == ontology) &
        (all_scores['output_type'] == 'RNA_SEQ')
    ]
    
    # Group by gene and compute mean absolute effect
    gene_impact = (
        filtered
        .groupby('gene_name')['raw_score']
        .apply(lambda x: x.abs().mean())
        .reset_index()
        .rename(columns={'raw_score': 'gene_impact'})
        .sort_values('gene_impact', ascending=False)
    )
    
    return gene_impact


def compute_rbi(gene_impact_df):
    """
    Compute patient-level Regulatory Burden Index.
    
    RBI = sum of all gene impacts / number of genes affected
    
    Normalization prevents bias toward patients with more variants.
    
    Args:
        gene_impact_df: DataFrame from compute_gene_impact() with 'gene_impact' column
    
    Returns:
        float: Normalized RBI score
    """
    # Sum all gene impacts
    total_burden = gene_impact_df['gene_impact'].sum()
    
    # Normalize by number of genes
    n_genes = len(gene_impact_df)
    
    # Return normalized RBI
    rbi = total_burden / n_genes if n_genes > 0 else 0.0
    
    return rbi
