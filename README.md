# Symmetric Non-negative Matrix Factorization (SymNMF) Clustering

This project was implemented as the final assignment in the Operating Systems / Software Project course. It includes a full implementation of the **SymNMF clustering algorithm** along with performance comparison to **K-Means** clustering.

## 📁 Contents

- \`symnmf.py\`: Python interface for all algorithm stages.
- \`symnmf.c\`, \`symnmf.h\`: C implementation of matrix computations.
- \`symnmfmodule.c\`: Python C API wrapper (for integration with Python).
- \`symnmf.py\`: Handles CLI arguments, matrix computation, and clustering via SymNMF.
- \`analysis.py\`: Performance evaluation and comparison between SymNMF and KMeans using silhouette score.
- \`setup.py\`: Builds the shared object (\`.so\`) file for Python module.
- \`Makefile\`: Compiles the C part into an executable.
- Additional utility C headers and source files as needed.

## 🔍 Project Overview

The goal was to implement **Symmetric Non-negative Matrix Factorization (SymNMF)** and use it for clustering data points. The implementation includes:

1. **Similarity matrix** construction.
2. **Diagonal degree matrix** and **normalized similarity matrix**.
3. **SymNMF optimization** using multiplicative updates.
4. Converting the result matrix \`H\` into discrete cluster assignments.
5. Integration of Python and C for performance and modularity.

The final deliverable supports multiple modes:
- \`sym\`: compute similarity matrix
- \`ddg\`: compute diagonal degree matrix
- \`norm\`: compute normalized similarity matrix
- \`symnmf\`: perform the full clustering pipeline and output matrix \`H\`

Example:

\`\`\`
python3 symnmf.py 3 symnmf input.txt
\`\`\`

## 📊 Analysis

The \`analysis.py\` file compares the quality of SymNMF against KMeans clustering using the \`silhouette_score\` from \`sklearn.metrics\`.

Example:

\`\`\`
python3 analysis.py 3 input.txt
\`\`\`

Output:
\`\`\`
nmf: 0.1162
kmeans: 0.1147
\`\`\`

## ⚙️ Build Instructions

To build the C extension:

\`\`\`
python3 setup.py build_ext --inplace
\`\`\`

To build the C executable:

\`\`\`
make
\`\`\`

## 🧾 Notes

- All SymNMF logic, C extension integration, and testing scripts were implemented by me.
- Testers for validation were also written by me.

## 📎 License

This repository is intended for **academic** and **portfolio** use only."
