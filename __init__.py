"""
OncoReg: Regulatory Burden Index for Cancer Patients

A Python package for computing regulatory burden in cancer using AlphaGenome.
"""

import os
import pandas as pd
from alphagenome.data import gene_annotation
from alphagenome.models import dna_client

__version__ = "0.1.0"

# Global state for AlphaGenome client and gene annotations
_dna_model = None
_gtf_transcript = None
_configured = False


def configure(api_key=None):
    """
    Configure OncoReg with your AlphaGenome API key.
    
    This must be called before using any scoring functions.
    
    Args:
        api_key (str, optional): Your AlphaGenome API key. 
                                 If not provided, will look for ALPHAGENOME_API_KEY 
                                 environment variable.
    
    Example:
        >>> import oncoreg
        >>> oncoreg.configure(api_key='your_key_here')
        >>> rbi = oncoreg.score_patient('patient.vcf', cancer_type='breast')
    
    Raises:
        ValueError: If no API key is provided or found in environment
    """
    global _dna_model, _gtf_transcript, _configured
    
    # Get API key from argument or environment
    if api_key is None:
        api_key = os.environ.get('ALPHAGENOME_API_KEY')
    
    if api_key is None:
        raise ValueError(
            "No API key provided. Either pass api_key argument or set "
            "ALPHAGENOME_API_KEY environment variable.\n\n"
            "Get your free API key at: https://deepmind.google.com/science/alphagenome"
        )
    
    print("Configuring OncoReg...")
    
    # Initialize AlphaGenome client
    print("  Connecting to AlphaGenome...")
    _dna_model = dna_client.create(api_key)
    print("  ✅ Connected to AlphaGenome")
    
    # Load gene annotations
    print("  Loading gene annotations (GENCODE v46, hg38)...")
    gtf = pd.read_feather(
        'https://storage.googleapis.com/alphagenome/reference/gencode/'
        'hg38/gencode.v46.annotation.gtf.gz.feather'
    )
    _gtf_transcript = gene_annotation.filter_protein_coding(gtf)
    _gtf_transcript = gene_annotation.filter_to_mane_select_transcript(_gtf_transcript)
    print(f"  ✅ Loaded {len(_gtf_transcript)} gene annotations")
    
    _configured = True
    print("\n🎉 OncoReg configured successfully!")


def _check_configured():
    """Internal function to ensure configure() has been called."""
    if not _configured:
        raise RuntimeError(
            "OncoReg has not been configured. Please call oncoreg.configure(api_key='your_key') first."
        )


def _get_dna_model():
    """Get the configured AlphaGenome model client."""
    _check_configured()
    return _dna_model


def _get_gtf():
    """Get the loaded gene annotation GTF."""
    _check_configured()
    return _gtf_transcript


# Import main user-facing functions from other modules
from .scoring import score_variant
from .aggregation import compute_gene_impact, compute_rbi
from .data import TISSUE_MAP

# High-level user API
def score_patient(variants, cancer_type='breast', return_details=False):
    """
    Compute Regulatory Burden Index for a single cancer patient.
    
    Args:
        variants: List of variants as tuples (chrom, pos, ref, alt)
                  OR path to VCF file (future support)
        cancer_type: Cancer type - 'breast', 'lung', 'colon', etc
        return_details: If True, return dict with RBI + gene impacts + metadata
    
    Returns:
        float: Patient RBI score
        OR dict: {'rbi': float, 'n_variants': int, 'n_genes': int, 'gene_impacts': DataFrame}
    
    Example:
        >>> import oncoreg
        >>> oncoreg.configure(api_key='your_key')
        >>> variants = [('chr17', 43044000, 'C', 'T'), ('chr17', 43050000, 'G', 'A')]
        >>> rbi = oncoreg.score_patient(variants, cancer_type='breast')
        >>> print(f"RBI: {rbi:.4f}")
    """
    _check_configured()
    
    # TODO: Add VCF parsing when isinstance(variants, str)
    if isinstance(variants, str):
        raise NotImplementedError("VCF file parsing not yet implemented. Pass list of variants.")
    
    # Score each variant
    print(f"Scoring {len(variants)} variants...")
    all_scores = []
    for i, (chrom, pos, ref, alt) in enumerate(variants):
        print(f"  Variant {i+1}/{len(variants)}: {chrom}:{pos} {ref}>{alt}")
        scores = score_variant(chrom, pos, ref, alt)
        all_scores.append(scores)
    
    # Aggregate into gene impacts
    gene_impact = compute_gene_impact(all_scores, cancer_type=cancer_type)
    
    # Compute RBI
    rbi = compute_rbi(gene_impact)
    
    # Return results
    if return_details:
        return {
            'rbi': rbi,
            'n_variants': len(variants),
            'n_genes': len(gene_impact),
            'gene_impacts': gene_impact,
            'total_burden': gene_impact['gene_impact'].sum()
        }
    else:
        return rbi


__all__ = [
    'configure',
    'score_patient',
    'score_variant',
    'compute_gene_impact',
    'compute_rbi',
    'TISSUE_MAP',
    '__version__',
]
