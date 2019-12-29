# DRjCC
Clustering scRNA-seq by joint dimension reduction and cell clustering


Overview: This is code to do joint learning dimension reduction and cell clustering of single-cell RNA-seq data given in the "experiment" section of the paper: Wenming Wu, Xiaoke Ma*. "Joint learning dimension reduction and cell clustering of single-cell RNA-seq data." The coding here is a generalization of the algorithm given in the paper

Files: drjcc.m - The main function.

main_drjcc.m - A script with a real scRNA-seq data to show how to run the code.

Data_zheng.mat - A real scRNA-seq data used in the cell type clustering example. We retain the genes that are expressed in at least 10 cells for the dataset. The data Zheng dataset contains 500 human peripheral blood mononuclear cells (PBMCs) sequenced using GemCode platform, which consists of three cell types, CD56+ natural killer cells, CD19+ B cells and CD4+/CD25+ regulatory T cells.

constructW.m - Compute adjacent matrix W.

NormalizeUV.m - Normalize data.

bestMap.m - permute labels of L2 to match L1 as good as possible.

compute_NMI.m - Program for calculating the Normalized Mutual Information (NMI) between two clusterings.

AMI.m - Program for calculating the Adjusted Mutual Information (AMI) between two clusterings.

ARI.m - Program for calculating the Adjusted Rand Index ( Hubert & Arabie) between two clusterings.

hungarian.m - Solve the Assignment problem using the Hungarian method.

data$Zhengexpr.csv - Zheng dataset original expression data. The original data can be downloaded from 10X GENOMICS website.

data$Zheng.celltype.csv - The cell type of Zheng dataset.

Contact:

Please send any questions or found bugs to Xiaoke Ma xkma@xidian.edu.cn
