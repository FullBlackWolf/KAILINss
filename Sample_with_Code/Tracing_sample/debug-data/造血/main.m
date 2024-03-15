clear
%% Parameters
% Crop between adjacent window edges
crop = 150;    
% Window moving step size
step = 200;
% Window size
windows_size = 500;
% And adjust the number of eigenvector in the sub-function eig_classify to 
% adjust the number of split and the refinement of split.

% The deduplication of window edges needs further improvement.


%% Load data  and  Split to compute
MM0 = load('D:\A Novel Neural Network by Topological Dynamic Space Transform Method\KAILIN\chen-其他数据（修正bladder）\chen-其他数据（修正bladder）\造血干细胞（只测prog）/similarity_matrix_full_blood.mat');
count_=readtable('D:\A Novel Neural Network by Topological Dynamic Space Transform Method\KAILIN\chen-其他数据（修正bladder）\chen-其他数据（修正bladder）\造血干细胞（只测prog）/similarity_matrix_full_blood.csv');
MM0 = MM0.similarity_matrix_full_blood;

%MM0=matrix;
MM0=MM0(contains(count_.label,'_prog'),contains(count_.label,'_prog'));
count_=count_(contains(count_.label,'_prog'),:);


%重新排序
[p,splitlist] = binary_corr_sorting(MM0);

MM=MM0(p,p);

split=[];



% %% Draw
% hold on
% imagesc(MM,[0.0,0.002]);colormap jet;axis equal
% for sp = 1:length(all_split)
%     sp = all_split(sp);
%     plot([sp,sp],[-10,length(MM)+10],'r')
%     plot([-10,length(MM)+10],[sp,sp],'r')
% end

% saveas(gcf, 'figure.pdf');

count_result=count_(p,:);
%writetable(count_result,'C:\Users\Administrator\Desktop\chenkailin2\造血干细胞（只测prog）\chensresult.csv')

%% 划分矩阵块 lyw
split_simple=[1 168 254 415 627 890 1123 1360 1578 1868 2260 2397 2556 2592 2690 2812 2960 3059 3185 3393];
simple_matrix=zeros(length(split_simple)-1,length(split_simple)-1);
for i=1:length(split_simple)-1
    for j=1:length(split_simple)-1
        simple_matrix(i,j)=mean(mean(MM(split_simple(i):split_simple(i+1)-1,split_simple(j):split_simple(j+1)-1)));
    end
end


simple_label = cell(length(split_simple)-1,1);
for i=1:length(split_simple)-1
    cell_array=count_result.label(split_simple(i):split_simple(i+1)-1);
    % 统计每个元素的出现次数
    [unique_elements, ~, element_indices] = unique(cell_array);
    element_counts = histcounts(element_indices, 'BinMethod', 'integers');
    
    % 找到出现次数最多的元素
    [max_count, max_count_index] = max(element_counts);
    most_frequent_element = unique_elements(max_count_index);

    if strcmp(cell2mat(most_frequent_element),'Monocyte_prog') ||strcmp(cell2mat(most_frequent_element),'Monocyte')
        simple_label(i) = {'Monocyte_bias'};
    else 
        simple_label(i) = {'neutrophil_bias'};
    end
end
%writecell(simple_label,'C:\Users\Administrator\Desktop\chenkailin2\造血干细胞（只测prog）\simple_label.csv')
%writematrix(simple_matrix,'C:\Users\Administrator\Desktop\chenkailin2\造血干细胞（只测prog）\simple_matrix.csv')



%--------------Methods 2------------
% 计算精度
correct =0;
all=0;
for i=1:length(split_simple)-1
    cell_array=count_result.label(split_simple(i):split_simple(i+1));
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

%{
M0 = MM;
for s =1:20
M0 = correlation_phi(M0);
end
M0(M0<0)=-1;M0(M0>0)=1;
[u,s,v]=svds(M0,3);
[~,p]=sort(u(:,1));
%figure(1)
subplot(2,3,1);imagesc(MM);
subplot(2,3,2);imagesc(MM(p,p));
subplot(2,3,3);imagesc(M0(p,p));
M1 = MM(p,p);
sp = sum(u(p,1)<0);
B1 = M1(1:sp,1:sp);
for s =1:20
B1 = correlation_phi(B1);
end
B1(B1<0)=-1;B1(B1>0)=1;
[u0,s0,v0]=svds(B1,3);
[~,p0]=sort(u0(:,1));
p1=p;
p1(1:sp)=p0;
B2 = M1(sp+1:end,sp+1:end);
for s =1:20
B2 = correlation_phi(B2);
end
B2(B2<0)=-1;B2(B2>0)=1;
[u0,s0,v0]=svds(B2,3);
[~,p0]=sort(u0(:,1));
p1(sp+1:end)=sp+p0;
%p1 = p0(p1);
%figure(2)
BM1 = blkdiag(B1,B2);
subplot(2,3,4);imagesc(MM);
subplot(2,3,5);imagesc(M1(p1,p1));
subplot(2,3,6);imagesc(BM1(p1,p1));


count_result=count_(p1,:);
writetable(count_result,'D:\A Novel Neural Network by Topological Dynamic Space Transform Method\测试三角函数卷积的相空间重构2023-8-15\测试三角函数卷积的相空间重构2023-8-15\final\chensresult_p1.xlsx')
%}

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
h.ColorLimits = [0, 0.0007]


%调换顺序
number_index=[12,13,4,1,7,8,9,2,3,     14,15,16,17,18,19,10,11,5,6,];
simple_label_result=simple_label_table(number_index);
simple_matrix_result=simple_matrix(number_index,number_index);
row_labels = simple_label_result;
column_labels = simple_label_result;
h = heatmap(simple_matrix_result);
h.YDisplayLabels = row_labels; % 设置行标签
h.XDisplayLabels = column_labels; % 设置列标签
h.ColorLimits = [0.0001, 0.0005]