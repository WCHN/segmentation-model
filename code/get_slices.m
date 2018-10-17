function [slices,miss,lblnDetbf] = get_slices(dat,model,opt,do,reg)

%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------

modality = dat.modality{1}.name; % Imaging modality

%--------------------------------------------------------------------------
% Get image data
%--------------------------------------------------------------------------

[obs,dm,~,~,scl] = get_obs(dat);

%--------------------------------------------------------------------------
% Missing data struct
%--------------------------------------------------------------------------

miss = get_par('missing_struct',obs);

%--------------------------------------------------------------------------
% Bias field
%--------------------------------------------------------------------------

if do.bf
    % Compute bias-field
    bf = get_bf(dat.bf.chan,dm);      
    
    lblnDetbf                   = bf;
    lblnDetbf(isnan(lblnDetbf)) = 1;        
    lblnDetbf                   = log(prod(lblnDetbf,2));  
    lblnDetbf                   = sum(lblnDetbf);
else  
    % Bias field not estimated
    bf = 1;    
    
    lblnDetbf = 0;
end

%--------------------------------------------------------------------------
% Labels (if provided)
%--------------------------------------------------------------------------

labels = get_labels(dat,opt); 

%--------------------------------------------------------------------------
% Warped template
%--------------------------------------------------------------------------

template = warp_template(model,reg.y,reg.Affine);  

%--------------------------------------------------------------------------
% Build struct holding slice data
%--------------------------------------------------------------------------

cl     = cell(dm(3),1);
slices = struct('obs',cl,'bf',cl,'template',cl,'code',cl,'bin_var',cl,'labels',cl);

for z=1:dm(3)
    slices(z) = gmm_img('getslice',z,dm,obs,bf,template,miss.C,labels,scl);
end
%==========================================================================