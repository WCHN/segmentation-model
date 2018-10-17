function Z = cluster2template(Z_cluster,part)
Iz  = size(Z_cluster,1);
lkp = part.lkp;
K   = max(lkp);
Z   = zeros([Iz K]);
for k=1:K
    Z(:,k) = sum(Z_cluster(:,lkp == k),2);
end
%==========================================================================