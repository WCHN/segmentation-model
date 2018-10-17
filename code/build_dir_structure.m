function [dat,opt] = build_dir_structure(dat,opt)

% Create universal 'write' directory
dir_output = opt.dir_output;
if exist(dir_output,'dir'), rmdir(dir_output,'s'); end; mkdir(dir_output);   
   
if opt.template.do
    % Create directory that will store all that is model related
    dir_model     = fullfile(dir_output,'model');
    mkdir(dir_model);
    opt.dir_model = dir_model;

    % Create directory that will store all initial velocities
    dir_vel     = fullfile(dir_output,'vel');
    mkdir(dir_vel);   
    opt.dir_vel = dir_vel;   

    % Create directory that will store template derivatives
    dir_a_der     = fullfile(dir_output,'a_der');
    mkdir(dir_a_der);
    opt.dir_a_der = dir_a_der;
    
    if opt.verbose.model >= 3
        % Create directory that will store 2D segmentations
        dir_seg2d     = fullfile(dir_output,'seg2d');
        mkdir(dir_seg2d);
        opt.dir_seg2d = dir_seg2d;
    end
else
    S0 = numel(dat);
    for s=1:S0
        nam = dat{s}.name;
        d   = fullfile(dir_output,nam);
        mkdir(d);      
        
        % Create directory that stores segmentations
        dir_seg0 = fullfile(d,'seg');
        mkdir(dir_seg0);         

        % Create directory that stores segmentations before cleaning up
        dir_seg_orig        = fullfile(dir_seg0,'orig');    
        mkdir(dir_seg_orig);      
        dat{s}.dir.seg_orig = dir_seg_orig;   

        % Create directory that stores final segmentations
        dir_seg        = fullfile(dir_seg0,'final');
        mkdir(dir_seg);
        dat{s}.dir.seg = dir_seg;                

        if opt.seg.write_mllabels
            % Create directory that ML estimates from segmentations
            dir_ml     = fullfile(dir_seg0,'ml');
            mkdir(dir_ml);
            dat{s}.dir.ml = dir_ml;   
        end

        % Create directory that stores initial velocities
        dir_def        = fullfile(d,'in-vel');
        mkdir(dir_def);
        dat{s}.dir.def = dir_def;   

        % Create directory that stores images
        dir_img        = fullfile(d,'img');
        mkdir(dir_img);
        dat{s}.dir.img = dir_img;   

        % Create directory that stores bias-field
        dir_bf        = fullfile(d,'bf');
        mkdir(dir_bf);
        dat{s}.dir.bf = dir_bf; 
        
        if S0 > 1
            % Create directory that stores temporary initial velocities
            dir_vel        = fullfile(d,'vel');
            mkdir(dir_vel);   
            dat{s}.dir.vel = dir_vel;   
        end
    end
end
%==========================================================================
