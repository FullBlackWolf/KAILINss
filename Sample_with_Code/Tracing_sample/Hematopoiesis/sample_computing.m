
function [simple_label,simple_matrix]=sample_computing(count_result,split_simple,MM)

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




end