%% Hematopoiesis
clc;clear;

%% Load data and Split to compute
MM0 = load('D:/A Novel Neural Network by Topological Dynamic Space Transform Method/KAILIN/Database/blood5k_matrix.mat');

%Input the Umap X and Y
count0_=readtable(['D:/A Novel Neural Network by Topological Dynamic Space Transform Method/KAILIN/Database/blood5k_matrix.csv']);

%Input the label
label0_=readtable(['D:/A Novel Neural Network by Topological Dynamic Space Transform Method/KAILIN/Database/blood5k_matrix_label.csv']);


MM0 = MM0.similarity_matrix_full_blood;

%compute the cluster
[p0,splitlist0] = binary_corr_sorting(MM0);

%Similarity matrix result
MM_result0 = MM0(p0,p0);

%Umap result
count_result0=count0_(p0,:);

%Label result
label_result0=label0_(p0,:);
final_result0=table();
final_result0(:,1) = count_result0(:,1);
final_result0(:,2) = count_result0(:,2);
final_result0(:,3) = label_result0(:,1);
split_result0=splitlist0;

%output data

%% Reprogramming-----------------------------------------------------------------------


%% Load data and Split to compute
MM1 = load('D:/A Novel Neural Network by Topological Dynamic Space Transform Method/KAILIN/Database/program2k_matrix.mat');

%Input the Umap X and Y
count1_=readtable(['D:/A Novel Neural Network by Topological Dynamic Space Transform Method/KAILIN/Database/program2k_matrix.csv']);

%Input the label
label1_=readtable(['D:/A Novel Neural Network by Topological Dynamic Space Transform Method/KAILIN/Database/program2k_matrix_label.csv']);


MM1 = MM1.similarity_matrix_full_sb_wsw;

%compute the cluster
[p1,splitlist1] = binary_corr_sorting(MM1);

%Similarity matrix result
MM_result1 = MM1(p1,p1);

%Umap result
count_result1=count1_(p1,:);

%Label result
label_result1=label1_(p1,:);
final_result1=table();
final_result1(:,1) = count_result1(:,1);
final_result1(:,2) = count_result1(:,2);
final_result1(:,3) = label_result1(:,1);
split_result1=splitlist1;

%output data


writetable(final_result0, 'D:/A Novel Neural Network by Topological Dynamic Space Transform Method/KAILIN/Database/blood5k_result.csv')
writetable(final_result1, 'D:/A Novel Neural Network by Topological Dynamic Space Transform Method/KAILIN/Database/program2k_result.csv')