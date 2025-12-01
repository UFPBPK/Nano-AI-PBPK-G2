# Nano-AI-PBPK-G2

# Second-Generation AI-Assisted PBPK Model for Nanoparticles 
### A second-generation artificial intelligence-assisted physiologically based pharmacokinetic model for predicting nanoparticle tumor delivery in mice

**Last Updated:** November 24, 2025  

---
## 📘 Overview

- This repository provides all curated datasets and model codes (PBPK modeling, machine learning, and the second-generation AI-assisted PBPK model) in the manuscript.
---

## 📂 Repository Structure

### **PBPK Modeling (PBPK)**
- EvaDat_1.csv: Observed nanoparticle concentration–time data.
- NanoPBPK.R: Main PBPK model code.
- Model calibration.R: Parameter estimation workflow using observed data.
- Data_PBPK parameters.csv: Final calibrated PBPK parameters.
---
### **Machine Learning (ML_AI)**
- Data-ML development.csv: Dataset used to train ML models for PBPK parameter prediction.
- Six individual ML parameter prediction code, such as Final- 2nd-AI_QSAR_KTRES50.ipynb.
---
### **AI-assisted PBPK (AI-PBPK)**
- data_total.csv: ML-predicted tumor-related parameters.
- EvaDat_1.csv: Observed nanoparticle concentration–time data.
- AI-assisted PBPK.R: End-to-end script for prediction and comparison with observed values.
---

## 🌟 Key Features

- **Gradient-based parameter search strategy** for rapid estimation of tumor-specific PBPK parameters.

- **Direct nanoparticle tumor delivery prediction** powered by AI-assisted PBPK model

- **Experimentally validated** using in-house animal data to assess predictive performance and model robustness.

- **Improved applicable scenarios over previous generation**, demonstrating stronger generalization and real-world applicability.
