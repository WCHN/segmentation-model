function show_segmentations(dat,opt)

figname = '(SPM) Segmentations';

% ---------------------------------------------------------------------
% Get figure (create if it does not exist)
f = findobj('Type', 'Figure', 'Name', figname);
if isempty(f)
    f = figure('Name', figname, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);  
clf(f);

populations = spm_json_manager('get_populations',dat);
P           = numel(populations);
S0          = numel(dat);

mx_rows = 16;
nrows   = min(S0,mx_rows); 
K       = opt.template.K;
ncols   = K + 2;

cnt_plots = 1;
for p=1:P
    population0 = populations{p}.name;
    modality    = populations{p}.type;
    
    cnt = 1;
    for s=1:S0
        population = dat{s}.population;

        if strcmp(population0,population)
            
            nii = nifti(dat{s}.pth.im2d);    
            img = single(nii.dat(:,:,:,:));
            
            if isfield(dat{s}.pth,'bfim2d')
                nii  = nifti(dat{s}.pth.bfim2d);    
                bfim = single(nii.dat(:,:,:,:));
                img  = [img', bfim'];
            else
                img  = img';
            end                        
            
            sb = subplot(nrows,ncols,[1:2] + (cnt_plots - 1)*ncols);
            if strcmpi(modality,'CT')
                imagesc(img,[0 100]); axis off xy;
            else
                imagesc(img); axis off xy;
            end
            colormap(sb,gray)
%             title(population)   

            nii = nifti(dat{s}.pth.seg2d);    
            Z   = single(nii.dat(:,:,:,:));
            
            img = [];
            for k=1:size(Z,4)
                img = [img Z(:,:,:,k)'];
            end

            sb = subplot(nrows,ncols,[3:ncols] + (cnt_plots - 1)*ncols);
            imagesc(img,[0 1]); axis off xy;             
            colormap(sb,pink)
        
            cnt_plots = cnt_plots + 1;
                        
            if cnt == ceil(mx_rows/P) || cnt_plots > mx_rows
                break
            end
            
            cnt = cnt + 1;
        end
    end
end

drawnow;
%==========================================================================