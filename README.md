# Profiling Blood Proteome in Autoimmune Disease

Olink PEA, 5 diseases

All functions are in the functions_panAI.R file. The analysis method is two-pronged, consisting of differential analysis (limma) and machine learning (glmnet lasso). which takes a list of diseases and returns a list of proteins that are associated with the diseases. The function uses Olink Explore PEA data. The data contains the proteins and the diseases they are associated with. Part of the larger Human Disease Blood Atlas project, as described in[ Álvez et al](https://www.science.org/doi/10.1126/science.adx2678). https://www.science.org/doi/10.1126/science.adx2678
