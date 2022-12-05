%==========================================================================
% HYCOM2FVCOM:
%   Select the nesting nodes and cells
%
% input  :
%   fgrid    --- the file containing the FVCOM grid (nc, grd, or 2dm)
%   fsigma   --- fvcom sigma file (dat)
%   fout     --- nesting grid output file name and path (mat)
%   obc_node --- fvcom open boundary node id
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
fgrid = './data/fvcom/Tokyo_bay_large_grd.dat';
fsigma = './data/fvcom/Tokyo_bay_large_sigma.dat';
fout = './output/Japan_nesting_grid.mat';
obc_node = 1 : 415;
n_layer = 2;
%--------------------------------------------------------------------------


% Read FVCOM grid and sigma
f = f_load_grid(fgrid);
sigma = read_sigma(fsigma);
f = f_calc_sigma(f, sigma);

% Find out the nesting boundary nodes and cells
[node_layer, cell_layer, node_weight, cell_weight] = f_find_nesting(f, obc_node, n_layer);
fn = f_load_grid_nesting(f, [node_layer{:}], [cell_layer{:}], 'Node_weight', [node_weight{:}], 'Cell_weight', [cell_weight{:}]);

% Save the grid
save(fout, 'fn');

% Draw figure
cm = cm_load('Blues', 'NColor', 5, 'Flip');
close all
figure
hold on
f_2d_range(f);
f_2d_coast(f);
f_2d_mesh(f, 'Color', [.7 .7 .7]);
plot(f.x(obc_node), f.y(obc_node), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 1.5)
% f_2d_cell(f, [cell_layer{:}])
f_2d_image(fn, fn.cell_weight);
caxis([0 1])
colormap(cm)
cb = colorbar('North');
cb.Position = [0.15 0.76 0.48 0.02]; 
xlabel('Longitude (^oE)')
ylabel('Latitude (^oN)')
mf_save('mesh_2_layers.png')
