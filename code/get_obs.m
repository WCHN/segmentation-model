function [obs,dm,mat,vs,scl,V,C,mn,mx,nam] = get_obs(varargin)

% Parse inputs
%--------------------------------------------------------------------------
p              = inputParser;
p.FunctionName = 'get_obs';
p.addRequired( 'dat',           @isstruct);
p.addParameter('do_scl', false, @islogical);
p.addParameter('mask',   true,  @islogical);
p.addParameter('val',    1000,  @isnumeric);
% p.addParameter('samp',   0,     @isnumeric);
p.parse(varargin{:});
dat            = p.Results.dat;
do_scl         = p.Results.do_scl;
val            = p.Results.val;
do_msk         = p.Results.mask;
% samp           = p.Results.samp;

% Sanity-check image data
%--------------------------------------------------------------------------
[dm,mat,vs,V,C] = check_obs(dat);

% % Sub-sampling    
% [subsmp,o,dm] = get_subsampling_grid(dm,vs,samp);

% Total number of voxels
I = prod(dm); 

modality = dat.modality{1}.name;
if isfield(dat.modality{1},'channel')
    % Multi-channel data
    %----------------------------------------------------------------------
    C   = numel(dat.modality{1}.channel);
    scl = zeros(1,C);
    obs = zeros(I,C,'single');
    for c=1:C
%         if samp > 0
%             obs1 = zeros(dm,'single');
%             for z=1:numel(subsmp.z0)
%                 obs1(:,:,z) = spm_sample_vol(V(c),subsmp.x0,subsmp.y0,subsmp.z0(z)*o,0);
%             end
%         else
            obs1 = dat.modality{1}.channel{c}.nii.dat(:,:,:);            
%         end
        
        if do_msk
            % Exclude locations where any of the images is not finite, or is zero. 
            msk        = spm_misc('msk_modality',obs1,modality);
            obs1(~msk) = NaN;
        end
        
        if strcmpi(modality,'ct')
            scl(c) = 1;
        else
            % Scaling factor to make intensities more similar
            scl(c) = double(val/nanmean(nanmean(nanmean(obs1(:,:,:)))));            
        end
        
        if do_scl
            obs(:,c) = scl(c)*obs1(:);
        else
            obs(:,c) = obs1(:);
        end
    end
    
    [~,nam] = fileparts(dat.modality{1}.channel{1}.nii.dat.fname);
else
    % Single-channel data
    %----------------------------------------------------------------------
%     if samp > 0
%         obs = zeros(dm,'single');
%         for z=1:numel(subsmp.z0)
%             obs(:,:,z) = spm_sample_vol(V,subsmp.x0,subsmp.y0,subsmp.z0(z)*o,0);
%         end
%     else
        obs = dat.modality{1}.nii.dat(:,:,:);
%     end
    
    if do_msk
        % Exclude locations where any of the images is not finite, or is zero. 
        msk       = spm_misc('msk_modality',obs,modality); % figure(666);imshow3D(msk)
        obs(~msk) = NaN;           
    end
    if strcmpi(modality,'ct')
        scl = 1;
    else
        % Scaling factor to make intensities more similar
        scl = double(val/nanmean(nanmean(nanmean(obs(:,:,:)))));
    end
    
    if do_scl
        obs = scl*single(obs(:));
    else
        obs = single(obs(:));
    end
    
    [~,nam] = fileparts(dat.modality{1}.nii.dat.fname);
end

mn = double(nanmin(obs,[],1));
mx = double(nanmax(obs,[],1));
%==========================================================================

%==========================================================================
function [dm,mat,vs,V,C] = check_obs(dat)

if isfield(dat.modality{1},'channel')
    V   = spm_vol(dat.modality{1}.channel{1}.nii.dat.fname);
    dm  = dat.modality{1}.channel{1}.nii.dat.dim;
    mat = dat.modality{1}.channel{1}.nii.mat;
    C   = numel(dat.modality{1}.channel);    
    for c=1:C
        dm1 = dat.modality{1}.channel{c}.nii.dat.dim;
        if ~isequal(dm,dm1)
           error('~isequal(dm,dm1)') 
        end
    end
else
    C   = 1;    
    V   = spm_vol(dat.modality{1}.nii.dat.fname);
    dm  = dat.modality{1}.nii.dat.dim;
    mat = dat.modality{1}.nii.mat;
end
vs = spm_misc('vxsize',mat);

if numel(dm)==2, dm(3) = 1; end

if any(dm == 1) && find(dm == 1) ~= 3
    error('find(dm == 1) ~= 3')
end
%==========================================================================