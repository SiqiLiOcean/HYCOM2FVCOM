%==========================================================================
% HYCOM2FVCOM:
%   Select the nesting nodes and cells
%
% input  :
%
%   fgrd     --- the FVCOM grd file (dat)
%   fdep     --- the FVCOM dep file (dat)
%   fsigma   --- the FVCOM sigma file (dat)
%   fgrid    --- the file containing the FVCOM grid and depth (nc, grd, or 2dm)
%    (use either fgrd+fdep+fsigma or fgrid containing both grid and depth)
%   fout     --- nesting grid output file name and path (mat)
%   obc_node --- the FVCOM open boundary node id
%   n_layer  --- number of nesting cell layers
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
addpath('~/tools/matFigure')

clc
clear

%--------------------------------------------------------------------------
% Input
fgrd = '../data/gom7_grd.dat';
fdep = '../data/gom7_dep.dat';
fsigma = '../data/gom7_sigma.dat';
fgrid = '';
fout = '../output/gom7_nesting_grid.mat';
obc_node = 1 : 130;
n_layer = 5;
%--------------------------------------------------------------------------


% Read FVCOM grid and sigma
if isempty(fgrid)
    [x, y, nv] = read_grd(fgrd);
    h = read_dep(fdep);
    f = f_load_grid(x, y, nv, h);
    sigma = read_sigma(fsigma);
    f = f_calc_sigma(f, sigma);
else
    f = f_load_grid(fgrid);
end

% Find out the nesting boundary nodes and cells
[node_layer, cell_layer, node_weight, cell_weight] = f_find_nesting(f, obc_node, n_layer);
fn = f_load_grid_nesting(f, [node_layer{:}], [cell_layer{:}], ...
                         'Node_weight', [node_weight{:}],     ...
                         'Cell_weight', [cell_weight{:}]);


% Save the grid
save(fout, 'fn');

% Draw figure
cm = cm_load('Blues', 'NColor', n_layer, 'Flip');
close all
figure
hold on
grid on
f_2d_range(f);
% f_2d_coast(f);
f_2d_mesh(f, 'Color', [.7 .7 .7]);
plot(f.x(obc_node), f.y(obc_node), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 1.5)
% f_2d_cell(f, [cell_layer{:}])
f_2d_image(fn, fn.cell_weight);
caxis([0.1 1.1])
colormap(cm)
cb = colorbar('North');
cb.Position = [0.40 0.30 0.48 0.02]; 
xlabel('X (1x10^2 km)')
ylabel('Y (1x10^2 km)')
mf_xtick_scale(gca, 1e-5)
mf_ytick_scale(gca, 1e-5)
title('Nesting node relaxation weight')
% mf_save('mesh_layers.png')


