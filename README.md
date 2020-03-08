# DRjCC
Clustering scRNA-seq by joint dimension reduction and clustering
Overview:

This is code to do joint learning dimension reduction and cell clustering of single-cell RNA-seq data given in the "experiment" section of the paper: 

Wenming Wu, Xiaoke Ma*. "Joint learning dimension reduction and cell clustering of single-cell RNA-seq data." 

The coding here is a generalization of the algorithm given in the paper. DRjCC is written in the MATLAB programming language. To use, please download the DRjCC folder and follow the instructions provided in the README.doc.

Files:

drjcc.m - The main function.

main_drjcc.m - A script with a real scRNA-seq data to show how to run the code.

Data_zheng.mat - A real scRNA-seq data used in the cell type clustering example.  We retain the genes that are expressed in at least 10 cells for the dataset. The data Zheng dataset contains 500 human peripheral blood mononuclear cells (PBMCs) sequenced using GemCode platform, which consists of three cell types, CD56+ natural killer cells, CD19+ B cells and CD4+/CD25+ regulatory T cells. 

constructW.m - Compute adjacent matrix W.

NormalizeUV.m - Normalize data.

bestMap.m - permute labels of L2 to match L1 as good as possible.

compute_NMI.m - Program for calculating the Normalized Mutual Information (NMI) between two clusterings.

AMI.m - Program for calculating the Adjusted Mutual Information (AMI) between two clusterings.

ARI.m - Program for calculating the Adjusted Rand Index ( Hubert & Arabie) between two clusterings.

hungarian.m - Solve the Assignment problem using the Hungarian method.

data$Zhengexpr.csv - Zheng dataset original expression data. The original data can be downloaded from 10X GENOMICS website.

data$Zheng.celltype.csv - The cell type of Zheng dataset. 

Example:

Follow the steps below to run DRjCC （also contained in the " main_drjcc.m" file）. Here use a real scRNA-seq data (Data_Zheng) set as an example.

clear all;
clc;
load('Data_zheng.mat')%Loading data
X=V;%%%%%%%Rows are genes, columns are cell sample
%==============Constructing a weight matrix==============
%Preset value before constructing weight matrix
options = [];
option.Metric = 'Cosine';
options.NeighborMode = 'KNN';%KNN
options.k =5;%5 nearest neighbors
options.WeightMode = 'Cosine';%Weights are 0 or 1, it can eplace with 'HeatKernel', 'Euclidean' 
 
W = constructW(X',options);
 
clear options;
options = [];
k1=100; k2=3;
%%%%%%%%%%%%%%k1 is the number of feature for dimension reduction,
%%%%%%%%%%%%%%k2 is the number of features for clustering in NMF(the number of cluster)
[Z_final,A_final,B_final, F_final] = drjcc(X, k1, k2, W, options); %% Call the main function to solve the variables
%%%%%%%%%%% Clustering cell type label
for e=1:size(F_final,2)
    v=F_final(:,e);
    ma=max(v);
    [s,t]=find(v==ma);
    l(e)=s;
end
%%%%%%%%%%%%%%=============Performance evaluation ================= 
ll=real_label;%%%  the label originally identified by the authors
l=l';%%% Labels obtained by DRjCC
[newl] = bestMap(ll,l);%% Permute label of l to match ll as good as possible
nmi=compute_NMI(ll,newl) %% Calculating the Normalized Mutual Information (NMI)
ami=AMI(ll,newl)%% Calculating the Adjusted Mutual Information (AMI)
ari = ARI(ll, max(ll),newl,max(newl)) %% Calculating the Adjusted Rand Index (ARI)

pre_label =newl;
if ~isempty(ll) 
    exact = find(pre_label == ll);
    accuracy = length(exact)/length(newl) %% Calculating the accuracy
    else
    accuracy = []
end


Contact:

Please send any questions or found bugs to Xiaoke Ma xkma@xidian.edu.cn 
