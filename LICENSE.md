# License and Data Usage Agreement

## DECODE Project - Academic Research Use Only

**Copyright © 2026 Prof. Dr. David Kaul**
**Charité – Universitätsmedizin Berlin**

---

## Overview

This repository contains analysis code and documentation for the **DECODE** project (DEconvolution of Epigenetic Composition in Oncological Diagnostics and Exploration). The project analyzes DNA methylation data from glioblastoma patients for academic research purposes.

---

## Code License

The analysis scripts, documentation, and software code in this repository are licensed under the **MIT License**:

```
MIT License

Copyright (c) 2026 David Kaul

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

## Data Restrictions

### ⚠️ CRITICAL: NO RAW DATA IS INCLUDED

**This repository does NOT contain any patient data or raw methylation files.**

The analysis pipeline is designed to work with:
1. **Public reference datasets** from Gene Expression Omnibus (GEO)
2. **Local Charité data** (NEVER committed to version control)

### Data Protection Compliance

All patient-derived methylation data is subject to:

- **GDPR** (General Data Protection Regulation)
- **Charité Data Protection Guidelines**
- **Ethics Committee Approval** (required before data access)

**Raw data files (.idat, .RData, .rds) are strictly prohibited from being committed to this repository.**

### Permitted Data

The following data types may be included in this repository:
- Aggregated, anonymized summary statistics
- Reference matrices derived from public GEO datasets
- Visualizations without patient identifiers
- Simulated/synthetic data for benchmarking

---

## Usage Restrictions

### Academic Use Only

This code is intended for:
- Academic research and education
- Non-commercial scientific analysis
- Dissertation work (HMU, Charité)

### Prohibited Uses

You may NOT use this code or derivatives for:
- Commercial diagnostic services
- Clinical decision-making without validation
- Any purpose violating patient privacy
- Redistribution of patient data

---

## Citation

If you use this code or methodology in your research, please cite:

```
Kaul, D. (2026). DECODE: Deconvolution of Epigenetic Composition in
Oncological Diagnostics and Exploration. Master's Thesis,
Charité – Universitätsmedizin Berlin.
```

---

## Third-Party Licenses

This project uses the following open-source software:

| Package | License | Source |
|---------|---------|--------|
| R | GPL-2 \| GPL-3 | https://www.r-project.org/ |
| Bioconductor | Artistic-2.0 | https://www.bioconductor.org/ |
| minfi | Artistic-2.0 | Bioconductor |
| deconvMe | GPL-3 | GitHub (omnideconv/methyldeconv) |
| ggplot2 | MIT | CRAN |

See individual package documentation for full license terms.

---

## Disclaimer

**NO WARRANTY**

This software is provided for research purposes only. It has NOT been validated for clinical use. The authors make no warranties about the accuracy, reliability, or suitability of this software for any particular purpose.

**RESEARCH USE ONLY - NOT FOR DIAGNOSTIC USE**

---

## Contact

**Principal Investigator:**
Prof. Dr. David Kaul
Email: dpk194@student.bham.ac.uk

**Supervisors:**
- Dr. Roland Arnold (Birmingham, UK)
- Prof. Dr. David Capper (Charité, Berlin)

---

## Acknowledgments

This project is part of a dissertation conducted at:
- **Charité – Universitätsmedizin Berlin**
- In collaboration with **University of Birmingham**

Funding and support details will be added upon project completion.

---

**Last Updated:** 2026-02-01
