function varargout = get_par(varargin)
%__________________________________________________________________________
% Get various parameters of the model
%
% FORMAT lkp     = get_par('lkp',modality,opt)
% FORMAT miss    = get_par('missing_struct',obs)
% FORMAT ix      = get_par('ix_zero_resp',population,lkp,opt)
% FORMAT tiny    = get_par('tiny')
% FORMAT [bg,gf] = get_par('bg_fg',bg,K)
%
% FORMAT help get_par>function
% Returns the help file of the selected function.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

if nargin == 0
    help get_par
    error('Not enough argument. Type ''help get_par'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch lower(id)
    case 'lkp'
        [varargout{1:nargout}] = get_lkp(varargin{:});      
    case 'missing_struct'
        [varargout{1:nargout}] = get_miss(varargin{:});          
    case 'ix_tiny'
        [varargout{1:nargout}] = get_ix_tiny(varargin{:});               
    case 'tiny'
        [varargout{1:nargout}] = get_tiny(varargin{:});          
    case 'bg_fg'
        [varargout{1:nargout}] = get_bg_fg(varargin{:});                
    otherwise
        help get_par
        error('Unknown function %s. Type ''help get_par'' for help.', id)
end
%==========================================================================

%==========================================================================
function [lkp,mg] = get_lkp(modality,opt)
% FORMAT [lkp,mg] = get_par('lkp',modality,opt)
%__________________________________________________________________________
K    = opt.template.K;
dict = opt.dict.lkp;
if dict.isKey(modality)
    lkp = dict(modality);    
else
    lkp = 1:K;
end

if max(lkp) ~= K    
    error('max(lkp) ~= K')
end

mg               = ones(1,numel(lkp));
for k=1:max(lkp)
    kk           = sum(lkp == k);
    mg(lkp == k) = 1/kk;
end   
%==========================================================================

%==========================================================================
function miss = get_miss(obs)
% FORMAT miss = get_par('missing_struct',obs)
%__________________________________________________________________________
miss    = struct;
miss.C  = spm_gmm_lib('obs2code',obs);
miss.L  = unique(miss.C);
miss.nL = numel(miss.L);
%==========================================================================

%==========================================================================
function [ix_tiny,ix0] = get_ix_tiny(population,lkp,opt)
% FORMAT ix_tiny = get_par('ix_tiny',population,lkp,opt)
%__________________________________________________________________________
K      = opt.template.K;
dict   = opt.dict.prop_excl;
if dict.isKey(population)
    ix0 = dict(population);    
else
    ix0 = zeros([1 K],'logical');
end

ix_tiny = zeros([1 numel(lkp)],'logical');
ix1     = find(ix0 == 1);

if ~isempty(ix1)
    for k=1:K
        if ismember(k,ix1)
           ix_tiny(lkp == k) = true;
        end
    end
end
%==========================================================================

%==========================================================================
function tiny = get_tiny
% FORMAT tiny = get_par('tiny')
%__________________________________________________________________________
tiny = 1e-4;
%==========================================================================

%==========================================================================
function [bg,fg] = get_bg_fg(bg,K)
% FORMAT [bg,fg] = get_par('bg_fg',bg,K)
%__________________________________________________________________________
fg  = 1:K;    
msk = ismember(fg,bg);
fg  = fg(~msk);
%==========================================================================