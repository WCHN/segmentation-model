function dat = update_prop(obs,bf,dat,template,labels,scl,dm,PropPrior,miss,opt)

% Parse input
cluster = dat.gmm.cluster;
armijo  = dat.armijo.prop;
part    = dat.gmm.part;
prop    = double(dat.gmm.prop(:))';
alpha   = double(PropPrior.alpha(:))';
% tiny    = get_par('tiny');

% Parameters
const    = spm_gmm_lib('Const', cluster{1}, cluster{2}, miss.L);
K        = size(template,2);
ix_tiny  = get_par('ix_tiny',dat.population,part.lkp,opt);

if numel(alpha) < numel(prop)
    alpha = padarray(alpha, [0 numel(prop) - numel(alpha)], 'post', 'replicate');
end
    
if armijo < 1e-6
    % Already found optimal solution
    dat.armijo.prop = min(armijo*1.25,1);
    
    return; 
end

%--------------------------------------------------------------------------
% GN loop
%--------------------------------------------------------------------------

for gn=1:opt.prop.niter 
    
    % Init
%     [dlb,mom] = gmm_img('init_lb_and_mom',miss);

    % Neighborhood part
    lnPzN = gmm_mrf('apply',dat.mrf);

    sumZ = 0; sumPi = 0; sumPi2 = 0;
    for z=1:dm(3)
        
        % Get slices
        [slice,ix] = gmm_img('getslice',z,dm,obs,bf,template,miss.C,labels,scl);

        if dat.mrf.do && numel(lnPzN) > K           
            lnPzNz = double(lnPzN(ix,:));
        else
            lnPzNz = lnPzN;
        end
    
        % Compute responsibilities and lb
%         [Z,dlb,BX] = gmm_img('slice_resp_and_lb',slice,cluster{1},cluster{2},prop,part,miss,const,lnPzNz,ix_tiny,dlb);
        Z = gmm_img('slice_resp_and_lb',slice,cluster{1},cluster{2},prop,part,miss,const,lnPzNz,ix_tiny);
        
%         % Compute sufficient statistics 
%         mom = gmm_img('slice_mom',mom,Z,slice,miss,BX);    
        
        % Map cluster responsibilities to tissue responsibilities
        Z = cluster2template(Z,part);   
        
        if dat.mrf.do        
            dat.mrf.oZ(:,:,z,:) = reshape(uint8((2^8)*Z),[dm(1:2) 1 K]);
        end
    
        [dsumZ, dsumPi, dsumPi2] = suffstat(prop, Z, slice.template);

        sumZ   = sumZ + dsumZ;
        sumPi  = sumPi + dsumPi;
        sumPi2 = sumPi2 + dsumPi2;
    end    

%     lbX = spm_gmm_lib('MarginalSum', mom.SS0, mom.SS1, mom.SS2, cluster{1}, cluster{2}, miss.L, mom.SS2b);    
% 
%     dat.lb.X(end + 1)       = lbX;     
%     dat.lb.Z(end + 1)       = dlb.Z; 
%     if numel(part.mg) > numel(prop)
%         dat.lb.mg(end + 1)  = dlb.mg;   
%     end
%     if ~isempty(labels)
%         dat.lb.lab(end + 1) = dlb.lab;
%     end
%     if dat.mrf.do   
%         dat.lb.ZN(end + 1)  = gmm_mrf('lowerbound',dat.mrf);
%     end
% 
%     dat.lb = check_convergence('gmm',dat.lb,1,opt.verbose.prop);
    
    p  = spm_matcomp('softmax',prop);    
    p2 = p'*p;

    % ---------------------------------------------------------------------
    % Compute grad/hess
    g = sumPi - sumZ - (alpha - 1).*(1 - p);
    H = diag(sumPi) - sumPi2 + diag(alpha - 1) .* (diag(p) - p2);
    H = spm_matcomp('LoadDiag', H);

    % ---------------------------------------------------------------------
    % GN step
    dw = H\g';
    dw = dw';

    % ---------------------------------------------------------------------
    % Line search    
    oprop  = prop;    
    for line_search=1:opt.nline_search.prop

        prop = oprop - armijo * dw;
        
        prop_reg = sum((alpha - 1) .* log(spm_matcomp('softmax',prop)));
        
        nlb = prop_reg;
        nlb = nlb + dat.lb.lab(end) + dat.lb.X(end) ...
              + dat.lb.MU(end) + dat.lb.A(end) + dat.lb.v_reg(end) ...
              + dat.lb.aff_reg(end) + sum(dat.lb.bf_reg(end,:)) + dat.lb.lnDetbf(end) + dat.lb.mg(end) + dat.lb.ZN(end);                        
        prop_reg = sum((alpha - 1) .* log(spm_matcomp('softmax',prop) + eps));

        [dlb,~,dat.mrf] = gmm_img('img_lb_and_mom',obs,bf,[],template,labels,prop,cluster{1},cluster{2},miss,part,dm,dat.mrf,ix_tiny,{'prop',oprop});                      
        nlb             = nlb + dlb.Z;             

        if nlb > dat.lb.last
            armijo = min(armijo*1.25,1);
            
            dat.lb.Z(end + 1)        = dlb.Z;  
            dat.lb.prop_reg(end + 1) = prop_reg;

            break;
        else
            armijo = armijo*0.5;

            if line_search == opt.nline_search.prop                        
                prop = oprop;            
            end        
        end

    end

    [dat.lb,gain] = check_convergence('prp',dat.lb,gn,opt.verbose.prop,armijo);
    
    if gain < opt.prop.tol
       % Finished updating registration
       break;
    end
end

if opt.verbose.prop >= 3
    % Plot GMM
    [MU,A] = get_mean_prec(cluster);
    spm_gmm_lib('Plot', 'GMM', obs, {MU,A}, spm_matcomp('softmax',prop), part);
end

dat.armijo.prop = armijo;
dat.gmm.prop    = prop;
%==========================================================================

%==========================================================================
function [sumZ, sumPi, sumPi2] = suffstat(w, Z, A)

K      = size(Z,2);
sumZ   = sum(double(Z),1);
Pi     = spm_matcomp('softmax',bsxfun(@plus, double(A), double(w)));
sumPi2 = zeros(K);
for k=1:K
    sumPi2(k,k) = sum(Pi(:,k).^2);
    for l = k+1:K
        sumPi2(k,l) = sum(Pi(:,k).*Pi(:,l));
        sumPi2(l,k) = sumPi2(k,l);
    end
end
sumPi = sum(Pi,1);
%==========================================================================