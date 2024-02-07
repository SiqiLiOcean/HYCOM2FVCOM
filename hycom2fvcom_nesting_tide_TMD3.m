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
% 2023-02-28  Siqi Li  Fixed the unit mismatch between TMD and HYCOM
%==========================================================================
addpath('~/tools/matFVCOM')
% addpath('~/tools/TMD')
addpath('~/tools/TMD3')
addpath('~/tools/t_tide')

clc
clear

%--------------------------------------------------------------------------
% Input
fnesting = '../output/ECS_nesting_grid_3d.mat';
% tide_name = ["M2" "N2" "S2" "K2" "K1" "O1" "P1" "Q1"];
tide_name = ["M2" "N2" "S2" "K2" "K1"];
Model = '/hosts/hydra.smast.umassd.edu/data5/siqili/data/TPXO9_atlas_v5/tpxo_binary2netcdf/TPXO9_atlas_v5.nc';
fout = '../output/ECS_nesting_tide.mat';
%--------------------------------------------------------------------------


% Load nesting grid
load(fnesting);

%--------------------------------------------------------------------------
% Extract the tide elevation components
% to check if all constituents are available
conList = tmd_conlist(Model);
for j = 1 : length(tide_name)
    k = find(ismember(upper(conList), tide_name{j}));
    if isempty(k)
        error(['Tide ' tide_name{j} ' is not included.'])
%     else
%         Cid(j) = k;
    end
end


%--------------------------------------------------------------------------
% Extract tide
for j = 1 : length(tide_name)
    amp(:,j) = squeeze(tmd_interp(Model, 'hAm', fn.y, fn.x, 'constituents', lower(tide_name{j}), 'coasts', 'unmask'));
    pha(:,j) = squeeze(tmd_interp(Model, 'hPh', fn.y, fn.x, 'constituents', lower(tide_name{j}), 'coasts', 'unmask')) / pi * 180;
end
% Conver phase from [-180 180] to [0 360]
pha = mod(pha, 360);
% Create the tidestruc for t_tide
for i = 1 : fn.node
    tide_zeta_struct(i,1) = create_tidestruc(tide_name, amp(i,:), pha(i,:));
end
clear amp pha

%--------------------------------------------------------------------------
% Extract the tide vector components
for j = 1 : length(tide_name)
    [fmaj(:,j), fmin(:,j), fpha(:,j), finc(:,j)] = tmd_ellipse(Model, lower(tide_name{j}), fn.yc, fn.xc);
end
% Create the tidestruc for t_tide
for i = 1 : fn.nele
    tide_uv_struct(i,1) = create_tidestruc(tide_name, fmaj(i,:), fmin(i,:), finc(i,:), fpha(i,:));
end
tide_uv_struct = repmat(tide_uv_struct, 1, fn.kbm1);


% Save out the file
save(fout, 'tide_zeta_struct', 'tide_uv_struct')

