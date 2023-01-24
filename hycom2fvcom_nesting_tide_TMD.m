%==========================================================================
% HYCOM2FVCOM:
%   Extract the tidal components for nesting boundary nodes and cells 
%   from OSU TPXO Tide Models with TMD
%
% input  :
%   fnesting --- fvcom nesting grid (mat, from hycom2fvcom_nesting_select.m)
%   tide_name--- tide constituent names (string array)
%   Model    --- tide model control file required by TMD, download from
%                https://www.esr.org/research/polar-tide-models/list-of-polar-tide-models/
%   fout     --- output tide coefficient file name and path (mat)
%
% output :
%   nesting grid mat file
%
% Siqi Li, SMAST
% 2022-11-30
%
% Updates:
%
%==========================================================================
addpath('~/tools/matFVCOM')
addpath('~/tools/TMD')
addpath('~/tools/t_tide')

clc
clear

%--------------------------------------------------------------------------
% Input
fnesting = './output/Japan_nesting_grid.mat';
tide_name = ["M2" "N2" "S2" "K2" "K1" "O1" "P1" "Q1"];
Model = '/hosts/hydra.smast.umassd.edu/data2/siqili/script/tools/TMD/DATA/Model_tpxo9';
fout = './output/Japan_nesting_tide.mat';
%--------------------------------------------------------------------------


% Load nesting grid
load(fnesting);

%--------------------------------------------------------------------------
% Extract the tide elevation components
[~, ~, ~, conList] = tmd_extract_HC(Model, fn.y(1), fn.x(1), 'z', []);
for j = 1 : length(tide_name)
    k = find(ismember(upper(conList), tide_name{j}, 'rows'));
    if isempty(k)
        error(['Tide ' tide_name{j} ' is not included.'])
    else
        Cid(j) = k;
    end
end

%--------------------------------------------------------------------------
% Extract tide
[amp, pha, depth] = tmd_extract_HC(Model, fn.y, fn.x, 'z', Cid);
amp = amp';
pha = pha';
% Create the tidestruc for t_tide
for i = 1 : fn.node
    tide_zeta_struct(i,1) = create_tidestruc(tide_name, amp(i,:), pha(i,:));
end
clear amp pha

%--------------------------------------------------------------------------
% Extract the tide vector components
for j = 1 : length(tide_name)
    [fmaj(:,j), fmin(:,j), pha(:,j), finc(:,j)] = tmd_ellipse(Model, fn.yc, fn.xc, tide_name{j});
end
% Create the tidestruc for t_tide
for i = 1 : fn.nele
    tide_uv_struct(i,1) = create_tidestruc(tide_name, fmaj(i,:), fmin(i,:), finc(i,:), pha(i,:));
end
tide_uv_struct = repmat(tide_uv_struct, 1, fn.kbm1);


% Save out the file
save(fout, 'tide_zeta_struct', 'tide_uv_struct')

