# macpie 1.1.0

## New functions
- Initial CRAN release of **macpie**: a scalable R toolkit for high-throughput transcriptomic (HTTr) analysis.
- Provides end-to-end support for MAC-seq and other HTTr data, including:
  - Quality control metrics:  `validate_metadata()`, `compute_qc_metrics()`, `plot_rle`, `plot_mds()`, etc.
  - Differential expression: `compute_single_de()`, `compute_multi_de()`, `plot_volcano()`, `plot_counts()`.
  - Support for chemical structure processing: `compute_chem_descriptors()`, `compute_smiles()`, `compute_single_dose_response()`, 
  `compute_multi_screen_profile()`.

