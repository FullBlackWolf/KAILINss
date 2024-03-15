clear


%% Load data and Split to compute
MM0 = load('D:/A Novel Neural Network by Topological Dynamic Space Transform Method/KAILIN/Database/program2k_matrix.mat');
%Input the Umap X and Y and label
count_=readtable(['D:/A Novel Neural Network by Topological Dynamic Space Transform Method/KAILIN/Database/重编程准备重排-matlab.csv']);
MM0 = MM0.similarity_matrix_full_sb_wsw;

%重新排序
[p,splitlist] = binary_corr_sorting(MM0);
MM=MM0(p,p);
count_result=count_(p,:);


%MM0=matrix;
% 选中prog。prog_index是_prog细胞的索引
MM_prog=MM(contains(count_.label,'_prog'),contains(count_.label,'_prog'));
%MM0=MM0(contains(count_(:,3),'_prog'),contains(count_(:,3),'_prog'));
count_prog=count_result(contains(count_result.label,'_prog'),:);





%% 划分矩阵块 lyw

split_simple=[1 196 667 827 1517 1741];
simple_matrix=zeros(length(split_simple)-1,length(split_simple)-1);
for i=1:length(split_simple)-1
    for j=1:length(split_simple)-1
        simple_matrix(i,j)=mean(mean(MM_prog(split_simple(i):split_simple(i+1)-1,split_simple(j):split_simple(j+1)-1)));
    end
end

% 统计每个划分块里面的
simple_label = cell(length(split_simple)-1,1);
for i=1:length(split_simple)-1
    cell_array=count_prog.label(split_simple(i):split_simple(i+1)-1);
    % 统计每个元素的出现次数
    [unique_elements, ~, element_indices] = unique(cell_array);
    element_counts = histcounts(element_indices, 'BinMethod', 'integers');
    
    % 找到出现次数最多的元素
    [max_count, max_count_index] = max(element_counts);
    most_frequent_element = unique_elements(max_count_index);

    if strcmp(cell2mat(most_frequent_element),'failed_prog') ||strcmp(cell2mat(most_frequent_element),'failed')
        simple_label(i) = {'failed_bias'};
    else 
        simple_label(i) = {'Reprogrammed_bias'};
    end
end
%writecell(simple_label,'C:\Users\Administrator\Desktop\chenkailin2\重编程（只测prog）\simple_label.csv')
%writematrix(simple_matrix,'C:\Users\Administrator\Desktop\chenkailin2\重编程（只测prog）\simple_matrix.csv')



%--------------Methods 2------------
% 计算精度
correct =0;
all=0;
for i=1:length(split_simple)-1
    cell_array=count_prog.label(split_simple(i):split_simple(i+1));
    cell_array=cell_array(contains(cell_array,'_prog'));
    % 统计每个元素的出现次数
    [unique_elements, ~, element_indices] = unique(cell_array);
    element_counts = histcounts(element_indices, 'BinMethod', 'integers');
    
    % 找到出现次数最多的元素
    [max_count, max_count_index] = max(element_counts);
    all = all+length(cell_array);
    correct = correct+max_count;
end

pred = correct/all
heatmap(simple_matrix)

simple_label_table=(simple_label);
% 创建行和列标签（示例）
row_labels = simple_label_table;
column_labels = simple_label_table;
% 使用 heatmap 函数并传递相应参数
h = heatmap(simple_matrix);
h.YDisplayLabels = row_labels; % 设置行标签
h.XDisplayLabels = column_labels; % 设置列标签




%调换顺序
number_index=[1,3,5,4,2];
simple_label_result=simple_label_table(number_index);
simple_matrix_result=simple_matrix(number_index,number_index);
row_labels = simple_label_result;
column_labels = simple_label_result;
h = heatmap(simple_matrix_result);
h.YDisplayLabels = row_labels; % 设置行标签
h.XDisplayLabels = column_labels; % 设置列标签