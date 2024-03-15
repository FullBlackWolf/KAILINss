
function [simple_label_result,simple_matrix_result]=cluster_map(simple_label,simple_matrix)

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
number_index=[13,2,3,6,8,7,1      ,12,9,10,11,4,5,15,16,17,14];
simple_label_result=simple_label_table(number_index);
simple_matrix_result=simple_matrix(number_index,number_index);
row_labels = simple_label_result;
column_labels = simple_label_result;
h = heatmap(simple_matrix_result);
h.YDisplayLabels = row_labels; % 设置行标签
h.XDisplayLabels = column_labels; % 设置列标签
h.ColorLimits = [0.00016, 0.0004]

end