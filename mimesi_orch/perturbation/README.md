# Emission Perturbation Module (MIMESI)

This module generates ensemble perturbations of emission fields for use in
DART-based chemical data assimilation experiments.

## Supported models
- FARM
- CHIMERE
- WRF-style NetCDF (curvilinear grids)

## Methodology
- Gaussian horizontal correlations
- Exponential vertical correlations
- AR(1) temporal correlation
- Ensemble mean preservation

## Configuration
All parameters are controlled via environment variables prefixed with:

EMISSION_

Example:
```bash
export EMISSION_MEMS=20
export EMISSION_CORR_LENGTH_HZ=100000
export EMISSION_SPREAD=1.6

