%==========================================================================
% HYCOM2FVCOM:
%   Create TS initial file from HYCOM datset
%
% input  :
%   fgrid    --- the file containing the FVCOM grid (nc, grd, or 2dm)
%   fhycom   --- hycom data (NetCDF)
%   fout     --- output initial TS file name and path (NetCDF)
% 
% output :
%   initial TS file
%
% Siqi Li, SMAST
% 2022-11-30
%
% Updates:
% 2023-01-03  Siqi Li  Convert the fvcom longitude into [0 360]
% 2023-01-04  Siqi Li  Used new interpolating method
%==========================================================================
addpath('~/tools/matFVCOM')

clc
clear

%--------------------------------------------------------------------------
% Input
fgrid = './data/fvcom/Tokyo_bay_large_grd.dat';
fhycom = './data/hycom/hycom_20020101_0000.nc';
fout = './output/Japan_ini_ts.nc';
%--------------------------------------------------------------------------


% Read FVCOM grid and sigma
f = f_load_grid(fgrid);

% Read HYCOM data
lon0 = ncread(fhycom, 'lon');
lat0 = ncread(fhycom, 'lat');
depth0 = ncread(fhycom, 'depth');
t0 = ncread(fhycom, 'water_temp');
s0 = ncread(fhycom, 'salinity');
time0 = ncread(fhycom, 'time')/24 + datenum(2000,1,1);

% Dimensions
nz0 = length(depth0);

% Interpolation
t = nan(f.node, nz0);
s = nan(f.node, nz0);
disp('Calculating interpolating weight.')
wh = interp_2d_calc_weight('GLOBAL_BI', lon0, lat0, f.x, f.y);
for iz = 1 : nz0
    disp(['Interpolating the ' num2str(iz) 'th layer of ' num2str(nz0) ' layers.'])
    % Horizontal interpolation
    t_layer = interp_2d_via_weight(t0(:,:,iz), wh);
    s_layer = interp_2d_via_weight(s0(:,:,iz), wh);
    % Fill the points that are recognize as land in HYCOM
    t_layer = f_fill_missing(f, t_layer);
    s_layer = f_fill_missing(f, s_layer);
    % Set the layer below the depth as nan
    k_bot = f.h <= depth0(iz);
    t_layer(k_bot) = nan;
    s_layer(k_bot) = nan;
    % Store the data
    t(:,iz) = t_layer;
    s(:,iz) = s_layer;
%     % t
%     t_layer = t0(:,:,iz);
%     kt = ~isnan(t_layer);
%     Ft = scatteredInterpolant(xx0(kt), yy0(kt), t_layer(kt), 'linear', 'nearest');
%     t(:,iz) = Ft(x360, f.y);
%     % s
%     s_layer = s0(:,:,iz);
%     ks = ~isnan(s_layer);
%     Fs = scatteredInterpolant(xx0(kt), yy0(kt), s_layer(kt), 'linear', 'nearest');
%     s(:,iz) = Fs(x360, f.y);
end

% Fill the points under the bathemetry
t = fillmissing(t, 'nearest', 2);
s = fillmissing(s, 'nearest', 2);

% % % % Plot (Uncomment this part for figures)
% % % figdir = '../output/';
% % % close all
% % % for iz = 1 : nz0
% % %     disp(['----' num2str(iz, '%2.2d')])
% % %     figure('Visible', 'off')
% % %     hold on
% % %     subplot(2,1,1)
% % %     f_2d_image(f, t(:,iz));
% % %     colorbar
% % %     caxis([0 35])
% % %    title([num2str(depth0(iz)) ' m (' num2str(iz, '%2.2d') '): Temperature'])
% % %     subplot(2,1,2)
% % %     f_2d_image(f, s(:,iz));
% % %     colorbar
% % %     caxis([25 40])
% % %     title([num2str(depth0(iz)) ' m (' num2str(iz, '%2.2d') '): Salinity'])
% % % 
% % %     ffig = [figdir '/init_ts_' datestr(time0, 'yyyymmddHH') '_layer_' num2str(iz,'%2.2d') '.png'];
% % %     mf_save(ffig);
% % %     close
% % % end


% Write initial TS output
write_initial_ts(fout, -depth0, t, s, time0);
