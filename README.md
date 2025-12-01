# Nano-AI-PBPK-G2

# Second-Generation AI-Assisted PBPK Model for Nanoparticles 
### A second-generation artificial intelligence-assisted physiologically based pharmacokinetic model for predicting nanoparticle tumor delivery in mice
The manuscript describing this model is under review. Please refer to the manuscript for details once it is published.

**Last Updated:** November 24, 2025  

---
## 📘 Overview

- This repository provides all curated datasets and model codes (PBPK modeling, machine learning, and the second-generation AI-assisted PBPK model) in the manuscript.
---

## 📂 Repository Structure

### **PBPK Modeling (PBPK)**
- Experiment.csv: Observed nanoparticle concentration–time data for the external validation.
- NanoPBPK.R: Main PBPK model code.
- Model calibration.R: Parameter estimation workflow using observed data.
- Data_PBPK parameters.csv: Final calibrated PBPK parameters.
- Linear regression: the code file to generate results in Figure 2.
---
### **Machine Learning (ML_AI)**
- Six individual ML parameter prediction code, such as Final- 2nd-AI_QSAR_KTRES50.ipynb.
---
### **AI-assisted PBPK (AI-PBPK)**
- AI-assisted PBPK.R: End-to-end script for prediction and comparison with observed values.
---

## 🌟 Key Features

- **Gradient-based parameter search strategy** for rapid estimation of tumor-specific PBPK parameters.

- **Direct nanoparticle tumor delivery prediction** powered by AI-assisted PBPK model

- **Experimentally validated** using in-house animal data to assess predictive performance and model robustness.

- **Improved applicable scenarios over previous generation**, demonstrating stronger generalization and real-world applicability.
