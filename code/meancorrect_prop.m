function dat = meancorrect_prop(dat,opt)
% FORMAT dat = meancorrect_prop(dat,opt)
% dat - Subjects data structure
% opt - Options structure
%
% Zero-centre proportion parameters across subjects.

if opt.prop.norm && opt.prop.do      
    S0 = numel(dat);
    
    p = 0;
    for s=1:S0
        p = p + dat{s}.gmm.prop;
    end
    
    p_avg = p/S0;
    clear r

    % Some verbose
    fprintf('p_avg = [');
    for i=1:numel(p_avg) - 1
        fprintf('%4.2f ',p_avg(i));
    end
    fprintf('%4.2f',p_avg(end))
    fprintf(']\n');
                
    for s=1:S0
        dat{s}.gmm.prop = dat{s}.gmm.prop - p_avg;
    end            
end
%==========================================================================