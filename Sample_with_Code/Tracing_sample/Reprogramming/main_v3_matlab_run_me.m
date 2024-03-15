clear;clc;
%% Parameters



%% Load data  and  Split to compute
%% Load data and Split to compute
MM0 = load('data/distance_matrix.mat');
MM0 = MM0.distance_matrix;

%MM0 = load('data/program2k_matrix.mat');
%MM0 = MM0.similarity_matrix_full_sb_wsw;

%Input the Umap XY and Label
count_=readtable(['data/merged_data.csv']);

%computing

MM0=MM0(contains(count_.label,'_prog'),contains(count_.label,'_prog'));
count_=count_(contains(count_.label,'_prog'),:);
[p,splitlist] = binary_corr_sorting(MM0);
[uniqueList, ~, ~] = unique(splitlist, 'stable');
MM=MM0(p,p);
split=[];
count_result=count_(p,:);
split_simple=uniqueList;
split_simple(1)=1;
split_simple=[split_simple,length(MM0)]
[simple_label,simple_matrix]=sample_computing(count_result,split_simple,MM);
[acc]=acc_computing(split_simple,count_result);
[predict_result]=predict_computing(count_,simple_label,split_simple);
count_result_out=[count_result,predict_result];
[simple_label_result,simple_matrix_result]=cluster_map(simple_label,simple_matrix);



writecell(simple_label_result,'result/Figure_c_label_blood.csv')
writematrix(simple_matrix_result,'result/Figure_c_matrix_blood.csv')
writetable(count_result_out,'result/umap_blood.csv')