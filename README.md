# OncoReg

**Regulatory Burden Index for Cancer Patients**

OncoReg computes a patient-level Regulatory Burden Index (RBI) by scoring non-coding somatic variants using AlphaGenome and aggregating their predicted effects on gene regulation.

## Installation

```bash
pip install oncoreg
```

## Quick Start

```python
import oncoreg

# Configure with your AlphaGenome API key
oncoreg.configure(api_key='your_alphagenome_key')

# Score a single patient
rbi = oncoreg.score_patient(
    vcf_path='patient_tumor.vcf',
    cancer_type='breast'
)

print(f"Patient RBI: {rbi:.4f}")
```

## Getting an AlphaGenome API Key

OncoReg requires access to the AlphaGenome API. Get your free API key at:
https://deepmind.google.com/science/alphagenome

The API is free for non-commercial research use.

## Supported Cancer Types

- `breast` - Breast cancer (uses breast epithelium tissue)
- `lung` - Lung cancer
- `colon` - Colon cancer

More cancer types will be added in future releases.

## What is RBI?

The Regulatory Burden Index (RBI) measures the cumulative disruption to gene regulation caused by non-coding somatic mutations in a tumor. 

Higher RBI indicates greater regulatory dysregulation, which may correlate with:
- Tumor aggressiveness
- Treatment response
- Patient survival outcomes

