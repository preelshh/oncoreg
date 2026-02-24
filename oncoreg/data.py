"""
Data constants and tissue mappings for OncoReg.
"""

# Map user-friendly cancer types to AlphaGenome ontology terms
TISSUE_MAP = {
    'breast': 'UBERON:0008367',   # breast epithelium
    'lung': 'UBERON:0002048',      # lung
    'colon': 'UBERON:0001157',     # colon
    'prostate': 'UBERON:0002367',  # prostate gland (placeholder - verify)
    'ovarian': 'UBERON:0000992',   # ovary (placeholder - verify)
}


def get_tissue_ontology(cancer_type):
    """
    Get AlphaGenome ontology term for a cancer type.
    
    Args:
        cancer_type: User-friendly cancer name ('breast', 'lung', etc)
    
    Returns:
        str: UBERON ontology term
    
    Raises:
        ValueError: If cancer type not supported
    """
    ontology = TISSUE_MAP.get(cancer_type.lower())
    if ontology is None:
        supported = ', '.join(TISSUE_MAP.keys())
        raise ValueError(
            f"Cancer type '{cancer_type}' not supported. "
            f"Supported types: {supported}"
        )
    return ontology
