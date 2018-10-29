function show_segmentations(dat)

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

mx_rows = 20;
nrows   = min(S0,mx_rows); 

cnt_plots = 1;
for p=1:P
    population0 = populations{p}.name;
        
    cnt = 1;
    for s=1:S0
        population = dat{s}.population;

        if strcmp(population0,population)

            nii = nifti(dat{s}.pth.seg2d);    
            Z   = single(nii.dat(:,:,:,:));
            
            img = [];
            for k=1:size(Z,4)
                img = [img Z(:,:,:,k)];
            end

            subplot(1,nrows,cnt_plots);
            imagesc(img',[0 1]); axis off xy;
            title(population)   
            
            cnt_plots = cnt_plots + 1;
                        
            if cnt == floor(mx_rows/P) || cnt_plots >= mx_rows
                break
            end
            
            cnt = cnt + 1;
        end
    end
end

drawnow;
%==========================================================================