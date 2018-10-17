function [MU,A] = get_mean_prec(cluster)
MU = cluster{1}{1};
K  = size(MU,2);
W0 = cluster{2}{1};
n0 = cluster{2}{2};
A  = bsxfun(@times,W0,reshape(n0,[1 1 K])); % E[Lambda]
%==========================================================================