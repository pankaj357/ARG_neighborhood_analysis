# ARG Neighborhood Analysis in ESKAPE Pathogens

This repository contains data and analysis scripts for the project:

**"Comparative Genomic Analysis of Antibiotic Resistance Gene Neighborhoods in ESKAPE Pathogens Reveals High-Risk Clusters and Mobile Genetic Element Associations"**

---

## 🧬 Project Overview

Antibiotic resistance among ESKAPE pathogens is a global public health threat. This study investigates:

- The clustering of antibiotic resistance genes (ARGs) within genomes.
- The proximity of ARGs to mobile genetic elements (MGEs), indicating potential for horizontal gene transfer.
- A comparative genomic analysis was performed across **12 strains** (2 per ESKAPE species) using **Prokka** and **AMRFinderPlus**, with custom Python scripts for neighborhood analysis.

---

## 📁 Repository Structure

```
ARG_neighborhood_analysis/
├── amrfinder_output/         # AMRFinderPlus .txt output files for all 12 strains
├── prokka_annotations/       # Prokka GFF3 annotation files for all 12 genomes
├── results/                  # Plots and summary .csv files
├── scripts/                  # Python scripts for ARG clustering and MGE proximity
├── supplementary/            # Supplementary table: strain name mapping
├── README.md                 # This file
```

---

## 🔬 Tools Used

- **Prokka v1.14.6** – for genome annotation  
- **AMRFinderPlus v4.0.23** – for antibiotic resistance gene detection  
- **Python 3.11** – for all custom analyses (arg clustering, MGE proximity)  
- **matplotlib / seaborn** – for figure generation  
- **pandas / Biopython** – for data parsing and computation  

---

## 📊 Key Figures

- `results/clustered_ARG_distribution.png`: Barplot showing clustered vs. non-clustered ARGs per strain
- `results/MGE_associated_ARGs.png`: Barplot showing MGE-associated vs. non-associated ARGs per strain

---

## 🧾 Supplementary Material

- `supplementary/Strain_Mapping_Table.csv`: Full strain names and corresponding figure labels (e.g., "Strain1", "Strain2")

---

## 📎 Citation

**Pankaj**  
*Indian Agricultural Research Institute (IARI), New Delhi*  
Manuscript submitted to peer-reviewed journal (under review).  
DO NOT cite until formally published.

---

## 📬 Contact

For questions or collaboration inquiries:

**Pankaj**  
Email: [ft.pank@gmail.com]  
Affiliation: Indian Agricultural Research Institute (IARI), New Delhi

---

## 📄 License

This repository is shared under the **MIT License** for academic and research use. See the `LICENSE` file for details.
