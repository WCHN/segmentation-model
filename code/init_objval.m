function [objval,dat] = init_objval(dat,P)
if nargin<2, P = 1; end

objval          = struct;
objval.main     = [];
objval.template = [];
for p=1:P
    objval.hyperpar(p).KL_qVpV = [];  
    objval.hyperpar(p).ElnDetV = [];
end    

S0  = numel(dat);
for s=1:S0
    [~,~,~,C,~,~,~,~] = obs_info(dat{s});

    lb0 = struct('sum', -Inf, 'last', -Inf, ...
                 'X', 0, 'lnDetbf', 0, 'Z', ...
                 0, 'MU', 0, ...
                 'A', 0, 'bf_reg', zeros(1,C), 'aff_reg', 0, 'v_reg', 0, 'lab', 0, 'prop_reg', 0, 'mg', 0, 'ZN', 0);

    dat{s}.lb = lb0;
end
%========================================================================== 