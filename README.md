# DECODE

**DEconvolution of Epigenetic Composition in Oncological Diagnostics and Exploration**

Reference-based deconvolution and in-silico purification of bulk DNA methylation
profiles, written for a thesis on the mesenchymal methylation class of
glioblastoma (MSc Bioinformatics, University of Birmingham, in cooperation with
the Institute of Neuropathology, Charite - Universitaetsmedizin Berlin).

A bulk methylation profile is an average over every cell in the sample. DECODE
estimates how much of it belongs to which cell type, and can then regress the
non-neoplastic contribution out of every CpG, leaving an estimate of the tumour
cells alone.

## What the pipeline does

1. **Reference panel.** 130 published profiles of sorted or cultured cell
   populations, covering 13 cell types: neurons, oligodendrocytes, astrocytes,
   microglia, six blood cell types, endothelial cells, fibroblasts and
   malignant GBM cells. Markers are selected one-versus-all (500 probes per cell
   type, 2,225 unique probes after deduplication).
2. **Deconvolution** with EpiDISH (RPC), giving a fraction per cell type.
3. **In-silico purification.** Per CpG, the observed value is regressed on four
   aggregated non-tumour fractions (blood-derived immune, microglia, neurons,
   glia); the tumour estimate is the prediction with all four set to zero.
4. **In-silico contamination.** The reverse experiment: a mesenchymal-like
   immune profile is mixed into immune-poor tumours in known proportions.

## Reference panel

| Cell type | Source | Samples |
|---|---|---|
| Blood, six types | FlowSorted.Blood.EPIC | 37 |
| Neurons | Hannon et al. | 28 |
| Oligodendrocytes | Hannon et al. | 24 |
| Astrocytes | Hannon et al. | 7 |
| Microglia | Hannon et al. | 17 |
| Endothelial | Moss et al. | 2 |
| Fibroblasts | dura fibroblasts | 10 |
| GBM malignant | Shen et al. (BTIC lines) | 5 |

130 profiles in total. Full provenance is in Table 2 of the manuscript.

## How well it works, and where it does not

| Test | Result |
|---|---|
| 13 control samples not used for reference construction | every cell type its own top prediction, mean 92.5% |
| 23 synthetic mixtures of known composition | r = 0.99, mean absolute error 1.2 percentage points |
| Tumour-free brain, 278 samples on 450k and 111 on EPIC | brain fraction 72 and 78%, malignant component 3.7 and 2.1% |

**Two limitations are load-bearing and are not hidden here.**

*Astrocytes and oligodendrocytes cannot be told apart in real tissue.* The
available sorted-nuclei references define astrocytes by exclusion, as nuclei
negative for NeuN, SOX10 and IRF8, so that fraction collects whatever is none of
the other three brain cell types. On tumour-free brain the two glial fractions
together match the published estimate (56% against 54%), but their split does
not. Glial fractions are therefore reported pooled throughout. Scripts 48 and 49
reproduce this.

*The absolute purity scale is uncertain.* DECODE returns the lowest tumour
fraction of the available estimates (median 31% against 54% for RF_Purify and
about 80% histologically). The fractions rank samples reliably; they are not
calibrated cell counts.

## Main result reproduced by this code

- Mesenchymal (MES) tumours carry the lowest tumour fraction (19.6%) and by far
  the highest immune content (33.3%, against 9.9 to 11.8% in the RTK classes).
- After purification of all 721,769 CpGs and reclassification with the Heidelberg
  classifier, **179 of 181 MES tumours lose the mesenchymal call** (98.9%): 169
  become RTK II, 5 become RTK I, 2 stay MES. The RTK classes are far more stable
  (RTK II 95%, RTK I 82%).
- Conversely, mixing an immune profile into ten immune-poor RTK tumours turns
  **all ten** into MES calls, dose-dependently, at an immune content that real
  MES tumours carry.

## Scope of this repository

This is the thesis subset of the DECODE working repository: the scripts that
reproduce the figures, tables and numbers of the dissertation, together with the
processed results they rest on. Exploratory work that is not part of the thesis
(MOSAIC single-cell integration, HLA methylation, the second contamination
round, the MGMT analyses, the TERT VAF comparison, the Berlin and Knoll cohorts)
is not included. The text of the dissertation is not published here either.

### Pipeline order

```
00 setup -> 01 reference panel -> 04/05 in-silico validation
10 cohort deconvolution -> 11/12 ITHresolveGBM projection
53 purification (4 components) -> [external classifier run] -> 14c call table
   13b/14b are the superseded 3-component run, kept as the sensitivity comparison
19/19b figures | 29 differentiation | 15b t-SNE | 58 manuscript numbers
31/32/33 in-silico contamination | 48/49 normal-brain validation
60 ITHresolveGBM profile-set sensitivity check
```

`scripts/_purification_config.R` decides in one place which purification the
analyses use (four components by default, three as the comparison).

### What the data folder contains, and what it does not

`data/processed/` holds the deconvolution fractions, classifier calls and derived
statistics behind the manuscript. Sample identifiers are replaced by running
numbers (S001 and upwards); the mapping is not published. Age, sex, free-text
diagnosis and all case, laboratory and sequencing identifiers are excluded, so
Table 1 of the manuscript cannot be rebuilt from this repository. Raw IDAT files
cannot be shared at all under GDPR and the Berlin Data Protection Act.

### Two steps cannot be run here

The Heidelberg Brain Tumor Classifier (mnp.v12epicv2) is not publicly available.
The purified and the contaminated MethylSet objects were classified at the
Institute of Neuropathology, Charite - Universitaetsmedizin Berlin, and only the
returned calls are included. Building the reference panel also needs the
published IDAT files of the sources listed above.

## Requirements

R 4.3 or newer, Bioconductor (EpiDISH, sesame, minfi, sva) and CRAN packages
(ggplot2, dplyr, tidyr, Rtsne, readxl). Exact versions are pinned in
`renv.lock`; `Rscript scripts/00_setup_environment.R --install` sets the
environment up.

## Author

David Kaul, MSc Bioinformatics, University of Birmingham.
Supervisors: Dr Roland Arnold (Birmingham), Prof. David Capper (Charite Berlin).

## License

See LICENSE.md

Generated by `tools/export_public.R` on 2026-09-29.
