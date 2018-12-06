function [dat,CM,opt] = build_label_cm(dat,opt)

% Parameters
S0 = numel(dat);

% For each subject with labels, build the label confusion matrix
for s=1:S0        
    
    CM = get_label_cm(dat{s}.population,opt);
    
    % Assign to dat
    dat{s}.gmm.cm = CM;
end
%===========================================================================