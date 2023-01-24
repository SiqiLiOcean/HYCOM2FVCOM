%==========================================================================
% HYCOM2FVCOM:
%   Write tide structure file of nesting nodes and cells
%
% input  :
%   fnesting   --- nesting grid structure (from hycom2fvcom_nesting_select.m)
%   ftide_zeta --- ASCII tidal zeta
%   ftide_uv   --- ASCII tidal uv
%   tide_name  --- tide name (string array)
%   tide_freq  --- tide frequency (number array)
%   fout       --- tide output file path and name
% 
% output :
%   tide output, containing tide_zeta_struct and tide_uv_struct
%
% Siqi Li, SMAST
% 2023-01-07
%
% Updates:
%
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
tide_name = ["M2" "N2" "S2" "K2" "K1" "O1" "P1" "Q1"];
tide_freq = [ 0.08051140; 0.07899925; 0.08333333; 0.08356149;
              0.04178075; 0.03873065; 0.04155259; 0.03721850];
fout = '../output/gom7_nesting_tide.mat';

% Read the fvcom nesting grid
load(fnesting);

% Creating tide name matrix
name = [];
for k = 1 : length(tide_name)
    name = [name; pad(tide_name{k}, 4)];
end

% Read the elevation components
tide_zeta = load(ftide_zeta);
tide_zeta = tide_zeta(:, 1:18);

% Creating the tide_zeta struct for t_tide
for i = 1 : fn.node
    tide_zeta_struct(i,1).name = name;
    tide_zeta_struct(i,1).freq = tide_freq;
    
    inode = fn.nesting_node(i);
    tide_zeta_struct(i,1).tidecon(:,1) = tide_zeta(inode, 3:2:end);
    tide_zeta_struct(i,1).tidecon(:,3) = tide_zeta(inode, 4:2:end);
    tide_zeta_struct(i,1).tidecon(:,2) = 0.002;
    tide_zeta_struct(i,1).tidecon(:,4) = 20;
end
    


% Read the vector componenets
tide_uv = load(ftide_uv);
tide_uv = tide_uv(:,7:end);
written_nele = size(tide_uv,1) / fn.kbm1;

% Creating the tide_uv struct for t_tide
for iz = 1:fn.kbm1
    for j = 1:fn.nele
        tide_uv_struct(j,iz).name = name;
        tide_uv_struct(j,iz).freq = tide_freq;

        tmp = tide_uv((iz-1)*written_nele+1:iz*written_nele,:);
        icell = fn.nesting_cell(j);
        tide_uv_struct(j,iz).tidecon(:,1) = tmp(icell,1:4:end);
        tide_uv_struct(j,iz).tidecon(:,2) = 0.002;
        tide_uv_struct(j,iz).tidecon(:,3) = tmp(icell,2:4:end);
        tide_uv_struct(j,iz).tidecon(:,4) = 0.002;
        tide_uv_struct(j,iz).tidecon(:,5) = tmp(icell,3:4:end);
        tide_uv_struct(j,iz).tidecon(:,6) = 20;
        tide_uv_struct(j,iz).tidecon(:,7) = tmp(icell,4:4:end);
        tide_uv_struct(j,iz).tidecon(:,8) = 20;
    end
end

save(fout, 'tide_zeta_struct', 'tide_uv_struct')


