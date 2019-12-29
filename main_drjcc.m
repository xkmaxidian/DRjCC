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
    %==============GNMF=================
    options = [];
    k1=500; k2=3;%%%%%%%%%%%%%%k1 is the number of feature for dimension reduction, k2 is the number of features for clustering in NMF
    [DT_final,V_final,B_final, F_final] = drjcc(X, k1, k2, W, options);
%%%%%%%%%%% Clustering cell type label
    for e=1:size(F_final,2)
    v=F_final(:,e);
    ma=max(v);
    [s,t]=find(v==ma);
    l(e)=s;
    end
    %%%%%%%%%%%%%%==================Performance evaluation===============================
    ll=real_label;
    l=l';
    [newl] = bestMap(ll,l);
    NMI=compute_NMI(ll,newl)
    [AMI_]=AMI(ll,newl)
    ARI = ARI(ll,k1,newl,k2)
    pre_label =newl;
    if ~isempty(ll) 
    exact = find(pre_label == ll);
    accuracy = length(exact)/length(newl)
    else
    accuracy = []
    end
    