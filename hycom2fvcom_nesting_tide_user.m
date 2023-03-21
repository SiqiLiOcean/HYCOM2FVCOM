%==========================================================================
% HYCOM2FVCOM:
%   Write tide structure file of nesting nodes and cells
%
% input  :
%   fnesting   --- nesting grid structure
%   ftide_zeta --- ASCII tidal zeta
%   ftide_uv   --- ASCII tidal uv
%   tides      --- tide name (string array)
%   fout       --- tide output file path and name
% 
% output :
%   tide output, containing tide_zeta_struct and tide_uv_struct
%
% Siqi Li, SMAST
% 2023-01-07
%
% Updates:
% 2023-02-02  Siqi Li  Use create_tidestruc to create tide structure
% 2023-03-21  Siqi Li  knnsearch -> ksearch
%==========================================================================
addpath('~/tools/matFVCOM')
addpath('~/tools/t_tide')

clc
clear

%--------------------------------------------------------------------------
% Input
fnesting = '../output/gom7_nesting_grid.mat';
ftide_zeta = '../data/tide_ellip_results/amppha_muti.dat';
ftide_uv = '../data/tide_ellip_results/tide_ellip/tide_uv.dat';
tides = ["S2" "M2" "N2" "K1" "O1" "K2" "P1" "Q1"];
fout = '../output/gom7_nesting_tide.mat';

% Read the fvcom nesting grid
load(fnesting);
% Get the longitudes and latitudes of nesting nodes and cells
% --- For Spherical Coordinate
% nesting_lon = fn.x;
% nesting_lat = fn.y;
% nesting_lonc = fn.xc;
% nesting_latc = fn.yc;
% --- For Cartisian Coordinate 
% (required to modified the following based on your own domain)
[nesting_lon, nesting_lat] = sp_proj('1802', 'inverse', fn.x, fn.y, 'm');
[nesting_lonc, nesting_latc] = sp_proj('1802', 'inverse', fn.xc, fn.yc, 'm');


% Read the elevation components
data = load(ftide_zeta);
x_zeta = data(:,1);
y_zeta = data(:,2);
[lon_zeta, lat_zeta] = sp_proj('1802', 'inverse', x_zeta, y_zeta, 'm');
tide_zeta = data(:, 3:18);

% Creating the tide_zeta struct for t_tide
id = ksearch([lon_zeta lat_zeta], [nesting_lon nesting_lat]);
for i = 1 : fn.node
%     inode = fn.nesting_node(i);
    inode = id(i);
    tide_zeta_struct(i,1) = create_tidestruc(tides, ...
                                             tide_zeta(inode,1:2:end), ...
                                             tide_zeta(inode,2:2:end));
end
    


% Read the vector componenets
data = load(ftide_uv);
lon_uv = data(:,3);
lat_uv = data(:,4);
tide_uv = data(:,7:end);
written_nele = size(data,1) / fn.kbm1;

% Creating the tide_uv struct for t_tide
idc = ksearch([lon_uv lat_uv], [nesting_lonc nesting_latc]);
for iz = 1:fn.kbm1
    disp(['---' num2str(iz)])
    for j = 1:fn.nele
%         icell = fn.nesting_cell(j);
        icell = idc(j);
        tmp = tide_uv((iz-1)*written_nele+1:iz*written_nele,:);
        tide_uv_struct(j,iz) = create_tidestruc(tides, ...
                                                tmp(icell,1:4:end), ...
                                                tmp(icell,2:4:end), ...
                                                tmp(icell,3:4:end), ...
                                                tmp(icell,4:4:end));
    end
end

save(fout, 'tide_zeta_struct', 'tide_uv_struct')


