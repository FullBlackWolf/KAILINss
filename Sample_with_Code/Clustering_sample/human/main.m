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
MM0 = load('data/matrix_human.mat');
count_=readtable('data/human.xlsx');
MM0 = MM0.matrix_human;
%MM0=matrix;
% 是否去除prog
% MM0=MM0(contains(count_.label,'_prog'),contains(count_.label,'_prog'));
% count_=count_(contains(count_.label,'_prog'),:);
%重新排序
[p,splitlist] = binary_corr_sorting(MM0);
[uniqueList, ~, ~] = unique(splitlist, 'stable');
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
%writetable(count_result,'chensresult.csv')


%% 划分矩阵块 lyw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%split_simple需要自己划
split_simple=uniqueList;
split_simple(1)=1;
split_simple=[split_simple,length(MM0)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

simple_matrix=zeros(length(split_simple)-1,length(split_simple)-1);
for i=1:length(split_simple)-1
    for j=1:length(split_simple)-1
        simple_matrix(i,j)=mean(mean(MM(split_simple(i):split_simple(i+1)-1,split_simple(j):split_simple(j+1)-1)));
    end
end

% 统计每个划分块里面的
simple_label = cell(length(split_simple)-1,1);
for i=1:length(split_simple)-1
    cell_array=count_result.label(split_simple(i):split_simple(i+1)-1);
    % 统计每个元素的出现次数
    [unique_elements, ~, element_indices] = unique(cell_array);
    element_counts = histcounts(element_indices, 'BinMethod', 'integers');
    
    % 找到出现次数最多的元素
    [max_count, max_count_index] = max(element_counts);
    most_frequent_element = unique_elements(max_count_index);

    % if strcmp(cell2mat(most_frequent_element),'failed_prog') ||strcmp(cell2mat(most_frequent_element),'failed')
    %     simple_label(i) = {'failed_bias'};
    % else 
    %     simple_label(i) = {'Reprogrammed_bias'};
    % end
    simple_label(i) = most_frequent_element;
end
%writecell(simple_label,'C:\Users\Administrator\Desktop\其他数据\mouse\simple_label.csv')
%writematrix(simple_matrix,'C:\Users\Administrator\Desktop\其他数据\mouse\simple_matrix.csv')




%--------------Methods 2------------
% 计算精度
correct =0;
all=0;
for i=1:length(split_simple)-1
    cell_array=count_result.label(split_simple(i):split_simple(i+1));
    % cell_array=cell_array(contains(cell_array,'_prog'));
    % 统计每个元素的出现次数
    [unique_elements, ~, element_indices] = unique(cell_array);
    element_counts = histcounts(element_indices, 'BinMethod', 'integers');
    
    % 找到出现次数最多的元素
    [max_count, max_count_index] = max(element_counts);
    all = all+length(cell_array);
    correct = correct+max_count;
end

pred = correct/all

%输出预测结果表格
predict_result=cell(height(count_),1);
simple_label=cell(simple_label);
for build_time = 1:1:length(split_simple)-1

    if build_time == 1

     predict_result(1:split_simple(2),1) = simple_label(build_time,1);


    else%if build_time == length(split_simple)-1

     predict_result(split_simple(build_time):split_simple(build_time+1))=simple_label(build_time,1);


    end

end


count_result_out=[count_result,predict_result];
writetable(count_result_out,'result/predict.csv')