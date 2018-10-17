function fct = get_facts(dat,model)

% Template
fct.templ.dm  = model.template.nii.dat.dim;   % Template orientation matrix   
fct.templ.mat = model.template.nii.mat;       % Template orientation matrix                  
fct.templ.K   = fct.templ.dm(4);

% Subject
[dm,mat,vs,C,nam,V,fnames,chn_names] = obs_info(dat);
modality                             = dat.modality{1}.name; 

chn_names{end + 1}        = 'Template';
chn_names{end + 1}        = 'Z';

fct.subj.dm  = dm;
fct.subj.mat = mat;
fct.subj.vs  = vs;
fct.subj.C   = C;
fct.subj.nam = chn_names;
fct.subj.I   = prod(dm(1:3));
fct.subj.mod = modality;
fct.subj.ff  = get_ff(vs); % Fudge factor (see spm_preproc8)   
%==========================================================================