
function [predict_result]=predict_computing(count_,simple_label,split_simple)

predict_result=cell(height(count_),1);
simple_label=cell(simple_label);
for build_time = 1:1:length(split_simple)-1

    if build_time == 1

     predict_result(1:split_simple(2),1) = simple_label(build_time,1);


    else%if build_time == length(split_simple)-1

     predict_result(split_simple(build_time):split_simple(build_time+1))=simple_label(build_time,1);


    end

end





end