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
MM0 = load('C:/Users/Administrator/Desktop/其他数据/human/matrix_human.mat');
count_=readtable('C:/Users/Administrator/Desktop/其他数据/human/human.xlsx');
MM0 = MM0.matrix_human;
%MM0=matrix;
% 是否去除prog
% MM0=MM0(contains(count_.label,'_prog'),contains(count_.label,'_prog'));
% count_=count_(contains(count_.label,'_prog'),:);
%重新排序
[p,splitlist] = binary_corr_sorting(MM0);

MM=MM0(p,p);

split=[];

for left=1:step:length(MM)
    right = min(left+windows_size,size(MM,2));
    sp=eig_classify(MM(left:right,left:right),false,false,false,3.5);
    if(left>1 && right<length(MM))
        sp = sp(((0+crop<sp).*(sp<windows_size-crop))==1);
    else
        if (left==1)
            sp = sp(sp<windows_size-crop);
        else
            sp = sp(0+crop<sp);
        end
    end
    split = [split,sp+left-1];
end


%% Concatenate and Deduplication
all_split=sort(split);
diff_all_split = diff(all_split);
left_split = [1,1+find(diff_all_split>1)];
right_split = [find(diff_all_split>1),length(all_split)];
for idx =1:length(left_split)
    all_split(idx)=mean(all_split(left_split(idx):right_split(idx)));
end
all_split(idx+1:end)=[];


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
writetable(count_result,'C:\Users\Administrator\Desktop\其他数据\human\chensresult.csv')

%% Draw
hold on
imagesc(MM,[0.0,0.002]);colormap jet;axis equal
for sp = 1:length(all_split)
    sp = all_split(sp);
    plot([sp,sp],[-10,length(MM)+10],'r')
    plot([-10,length(MM)+10],[sp,sp],'r')
end
%% 划分矩阵块 lyw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%split_simple需要自己划
split_simple=[1 79 95 137 157 175 306 464 652 914 1140 1356 1509 1667 1860 1938 1958 1985 2000 2050 2283 2505 2777 3025 3316 3524 3757 3979 4201 4471 4750 4764 5032 5257 5287 5430 5497 5532 5769 5952 6140 6352 6487 6688 6963 7054 7205 7355 7510 7892 8136 8345 8568];
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
writecell(simple_label,'C:\Users\Administrator\Desktop\其他数据\human\simple_label.csv')
writematrix(simple_matrix,'C:\Users\Administrator\Desktop\其他数据\human\simple_matrix.csv')



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