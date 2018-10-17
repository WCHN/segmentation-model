function set_globals

% Set random seed, in case random numbers are used
rng('default');
rng(1);

% Set boundary conditions (BOUND_CIRCULANT 0, BOUND_NEUMANN 1)
spm_diffeo('boundary',0); 
spm_field('boundary',1);
%==========================================================================