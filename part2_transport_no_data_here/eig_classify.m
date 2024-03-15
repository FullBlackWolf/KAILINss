% MM = load('C:/Users/czwy/Desktop/1091160905/FileRecv/matrix_blood.mat');
% MM = MM.matrix_zaoxue;
% Matrix = MM;
% 
% is_draw_matrix_and_split = true;
% is_draw_eigenvectors = false;
% is_draw_split_lines = true;
function [all_split]=eig_classify(Matrix,is_draw_matrix_and_split,is_draw_eigenvectors,is_draw_split_lines,deta)
%%   Eigenvector decomposition
%
% for left = 2001:500:size(Matrix,1)
% right = min(left+500,size(Matrix,2));
% sub_matrix = Matrix(left:right,left:right);
sub_matrix = Matrix;%(2001:2500,2001:2500);
Adj = gather(sub_matrix);
L = -spdiags(sum(Adj)',0,-Adj);
D = spdiags(-L,0);
D = diag(mean(D)+D);
L = D^(-1/2)*L*D^(-1/2);
[e,v]  = eigs(sparse(L),64,'sm', struct('disp',0));
% max_e = max(real(e),[],1);
% min_e = min(real(e),[],1);
% std_e = std(real(e),[],1);
%mean(e(:,1:4))
% break
% end

%%   Find eigenvectors and Filter maximum and minimum values
%
if(is_draw_matrix_and_split)
    scrsz = get(groot,'ScreenSize');
    figure('Position',[scrsz(3)*0.25 scrsz(4)*0.125 scrsz(3)*0.5 scrsz(4)*0.75]);
    hold on
    imagesc(sub_matrix);colormap jet;axis equal
end
g=0;
all_split = [];
for s =1:length(v)
    if(isreal(e(:,s)))
        g=g+1;
        
        if(is_draw_matrix_and_split && is_draw_eigenvectors)
            % draw eigenvectors
            plot(e(:,s)*length(sub_matrix)/10-length(sub_matrix)/25*g)
        end
        idx1 = abs(e(:,s)-mean(e(:,s)))>deta*std(e(:,s));
        se = e(idx1,s);
        idx2 = find(diff(diff(abs(se)))<0)+1;
        idx3 = find(idx1);
        idx4 = idx3(idx2);
        if ~isempty(idx4)
            for sp = 1:length(idx4)
                sp = idx4(sp);
                all_split = [all_split,sp];
                
                if(is_draw_matrix_and_split && is_draw_eigenvectors)
                    % draw maximum and minimum values on eigenvectors
                    plot(sp,e(sp,s)*length(sub_matrix)/10-length(sub_matrix)/25*g,'*')
                end
            end
        end
        %[s,max_e(s),min_e(s),std_e(s),v(s,s),sum(e(:,s))]
    end
end


%%  Deduplication, wiht two implementations
%
all_split = unique(sort(all_split));
%disp(all_split)
if(0) %
    for sp = 1:length(all_split)
        for i=sp:length(all_split)
            if i+1<=length(all_split)
                % disp([sp,i,all_split(i+1),all_split(i)])
                if all_split(i+1)>all_split(i)+1
                    % disp([sp,i,all_split(i+1),all_split(i)])
                    break
                else
                    continue
                end
                break
            end
        end
        all_split(sp:i)=mean(all_split(sp:i));
    end
    all_split = unique(sort(all_split));
else
    diff_all_split = diff(all_split);
    left_split = [1,1+find(diff_all_split>1)];
    right_split = [find(diff_all_split>1),length(all_split)];
    %     disp(all_split);
    %     disp(diff_all_split);
    %     disp(left_split);
    %     disp(right_split);
    for idx =1:length(left_split)
        all_split(idx)=mean(all_split(left_split(idx):right_split(idx)));
    end
    all_split(idx+1:end)=[];
end

%%  Draw split lines
%
%disp(all_split)
if(is_draw_matrix_and_split && is_draw_split_lines)
    for sp = 1:length(all_split)
        sp = all_split(sp);
        plot([sp,sp],[-10,length(sub_matrix)+10],'r')
        plot([-10,length(sub_matrix)+10],[sp,sp],'r')
    end
end

end