function [dat,bf] = update_bf(dat,obs,bf,template,dm,labels,miss,opt)
% FORMAT [dat,bf] = update_bf(dat,obs,bf,template,dm,labels,miss,opt)
% dat      - Subject's data structure (one subject)
% obs      - Observed image
% bf       - Bias fields
% template - Warped + softmaxed template
% dm       - Image dimensions
% labels   - Manual labels
% miss     - Missing data structure
% opt      - Options structure
%
% Update bias fields, one at a time, by Gauss-Newton optimisation.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% Parse input
chan    = dat.bf.chan;      % Bias field encoding (C,B1,B2,B3,T)
cluster = dat.gmm.cluster;  % GMM posterior parameters
prop    = dat.gmm.prop;     % Tissue proportions
armijo  = dat.armijo.bf;    % Previous armijo factor
verbose = opt.verbose.bf;   % Verbosity level
part    = dat.gmm.part;     % Clusters to template mapping (lkp,mg)

% Parameters
kron   = @(a,b) spm_krutil(a,b); % Redefine kronecker product
[MU,A] = get_mean_prec(cluster); % Get means and precisions
K      = size(MU,2);             % Number of clusters
C      = numel(chan);            % Number of chqnnels
const  = spm_gmm_lib('Const', cluster{1}, cluster{2}, miss.L); % GMM likelihood norm
ix_tiny = get_par('ix_tiny',dat.population,part.lkp,opt);      % Labels to template mapping

for c=1:C % Loop over channels

    if armijo(c) < 1e-6
        % Already found optimal solution
        armijo(c) = min(armijo(c)*1.25,1);
        dat.lb    = check_convergence('bf',dat.lb,c,verbose,armijo(c)); 
        
        continue; 
    end
    
    % Init
    [dlb,mom] = gmm_img('init_lb_and_mom',miss);
    
    % Neighborhood part
    lnPzN = gmm_mrf('apply',dat.mrf);

    %----------------------------------------------------------------------
    % Compute objective function and its first and second derivatives
    %----------------------------------------------------------------------
    
    d3 = numel(chan(c).T); % Number of DCT parameters
    H  = zeros(d3,d3); % Second derivatives
    gr = zeros(d3,1);  % First derivatives    
    
    for z=1:dm(3) % Loop over slices
        
        % Get slices
        [slice,ix] = gmm_img('getslice',z,dm,obs,bf,template,miss.C,labels);

        if dat.mrf.do && numel(lnPzN) > K           
            lnPzNz = double(lnPzN(ix,:));
        else
            lnPzNz = lnPzN;
        end
    
        % Compute responsibilities and lb
        [Z,dlb,BX] = gmm_img('slice_resp_and_lb',slice,cluster{1},cluster{2},prop,part,miss,const,lnPzNz,ix_tiny,dlb);
        
        if dat.mrf.do        
            dat.mrf.oZ(:,:,z,:) = reshape(uint8((2^8)*cluster2template(Z,part)),[dm(1:2) 1 max(part.lkp)]);
        end
        
        % Compute sufficient statistics 
        mom = gmm_img('slice_mom',mom,Z,slice,miss,BX);
        
        % Get mask 
        msk = isfinite(sum(slice.obs,2));
                
        % Check to see if there are observations in this slice
        nm = nnz(msk);
        if nm == 0, continue; end            
        
        % Get observations
        cr = cell(C,1);
        for c1=1:C % Loop over channels
            cr{c1} = BX(msk,c1); 
        end

        % Compute derivatives
        w1 = zeros(nm,1);
        w2 = zeros(nm,1);
        for k=1:K % Loop over classes                
            w0  = zeros(nm,1);
            for c1=1:C % Loop over channels
                w0 = w0 + A(c1,c,k)*(MU(c1,k) - cr{c1});
            end
            w1  = w1 + Z(msk,k).*w0;
            w2  = w2 + Z(msk,k)*A(c,c,k);
        end
        wt1       = zeros(dm(1:2));
        wt1(msk) = -(1 + cr{c}.*w1); % US eq. 34 (gradient)
        wt2       = zeros(dm(1:2));
        wt2(msk) = cr{c}.*cr{c}.*w2 + 1; % Simplified Hessian of US eq. 34
        clear cr

        b3 = chan(c).B3(z,:)';
        gr = gr + kron(b3,spm_krutil(wt1,chan(c).B1,chan(c).B2,0));
        H  = H  + kron(b3*b3',spm_krutil(wt2,chan(c).B1,chan(c).B2,1));
        clear wt1 wt2 b3 Z msk

    end % <-- Loop over slices

    lbX = spm_gmm_lib('MarginalSum', mom.SS0, mom.SS1, mom.SS2, cluster{1}, cluster{2}, miss.L, mom.SS2b);    

    dat.lb.X(end + 1)       = lbX;     
    dat.lb.Z(end + 1)       = dlb.Z; 
    if numel(part.mg) > numel(prop)
        dat.lb.mg(end + 1)  = dlb.mg;   
    end
    if ~isempty(labels)
        dat.lb.lab(end + 1) = dlb.lab;
    end
    if dat.mrf.do   
        dat.lb.ZN(end + 1)  = gmm_mrf('lowerbound',dat.mrf);
    end

    dat.lb = check_convergence('bf',dat.lb,1,verbose);
    
    % Inverse covariance of priors
    ICO = chan(c).C;         

    % Gauss-Newton update of bias field parameters
    Update = reshape((H + ICO)\(gr + ICO*chan(c).T(:)),size(chan(c).T));
    clear H gr

    % Line-search
    %------------------------------------------------------------------
    obf          = bf;
    oT           = chan(c).T;
    bf_reg       = dat.lb.bf_reg(end,:);
    nline_search = opt.nline_search.bf;
    for line_search=1:nline_search

        % Update bias-field parameters
        chan(c).T = chan(c).T - armijo(c)*Update;

        % Compute new bias-field (only for channel c)
        [bf,bf_reg] = get_bf(chan,dm,bf,c,bf_reg);                                                  

        % Compute new lower bound
        nlb = sum(bf_reg);        
        nlb = nlb + dat.lb.Z(end) + dat.lb.lab(end) ...
              + dat.lb.MU(end) + dat.lb.A(end) + dat.lb.v_reg(end) ...
              + dat.lb.aff_reg(end) + dat.lb.prop_reg(end) + dat.lb.mg(end) + dat.lb.ZN(end);                        

        lblnDetbf                   = bf;
        lblnDetbf(isnan(lblnDetbf)) = 1;        
        lblnDetbf                   = log(prod(lblnDetbf,2));  
        lblnDetbf                   = sum(lblnDetbf);
        nlb                         = nlb + lblnDetbf;

        [dlb,~,dat.mrf] = gmm_img('img_lb_and_mom',obs,bf,[],template,labels,prop,cluster{1},cluster{2},miss,part,dm,dat.mrf,ix_tiny,{'bf',obf});                      
        nlb             = nlb + dlb.X;

        % Check new lower bound
        if nlb > dat.lb.last               
            armijo(c) = min(armijo(c)*1.25,1);

            dat.lb.X(end + 1)        = dlb.X;                    
            dat.lb.lnDetbf(end + 1)  = lblnDetbf;  
            dat.lb.bf_reg(end + 1,:) = bf_reg;                                                                                                            
            
            break;
        else                                
            armijo(c) = armijo(c)*0.5;
            chan(c).T = oT;

            if line_search == nline_search                                        
                bf = obf;                                         
            end
        end
    end
    clear oT Update obf
    
    dat.lb = check_convergence('bf',dat.lb,c,verbose,armijo(c));
end % <-- Loop over channels

if verbose >= 3
    show_bf_and_ivel(obs,dm,bf);
end
            
% Get DC component
dc.ln = zeros(1,C);
for c=1:C
    dc.ln(c) = chan(c).T(1,1,1);
end
dc.int = scl_from_bf(chan);

% Set output
dat.bf.chan   = chan;
dat.bf.dc     = dc;
dat.armijo.bf = armijo;
%==========================================================================