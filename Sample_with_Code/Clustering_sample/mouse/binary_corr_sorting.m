function [p,splitlist] = binary_corr_sorting(M)
M0 = M;
for s =1:20
    M0 = correlation_phi(M0);
end
M0(M0<0)=-1;
M0(M0>0)=1;
%[u,s,~]=svds(M0,2);
u=M0(:,1:2);
[~,p]=sort(u(:,1));
split=sum(u(p,1)<0);

%% 
if length(M)<150 || split<5 || length(M)-split<5
    %disp([length(M),split])
    p=1:length(M);
    split=0;
    splitlist = [split];
else
    %disp([length(M),split])
    [p1,s1] = binary_corr_sorting(M(p(1:split),p(1:split)));
    [p2,s2] = binary_corr_sorting(M(p(split+1:end),p(split+1:end)));
    p = [p(p1);p(p2+split)];
    splitlist = [s1,split,split+s2];
end

    