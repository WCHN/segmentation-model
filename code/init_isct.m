function dat = init_isct(dat)
S0 = numel(dat);   
for s=1:S0
    if isfield(dat{s}.modality{1},'channel')
        C    = numel(dat{s}.modality{1}.channel);
        isct = false(1,C);
        for c=1:C
            name = dat{s}.modality{1}.channel{c}.name;
            if strcmpi(name,'CT')
                isct(c) = true;
            end
        end
    else
        isct = true;
    end
    
    dat{s}.isct = isct;
end
%==========================================================================