function pD=correlation_phi(D)
dik=D - mean(D,2)*ones(1,length(D));
Dij = dik*dik';
dij = sparse(diag(diag(Dij).^(-1/2)));
pD = dij*Dij*dij;