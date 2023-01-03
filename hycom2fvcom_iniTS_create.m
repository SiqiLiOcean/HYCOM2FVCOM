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
% Convert longitude from 180 to 360 (the same as HYCOM).
x360 = calc_lon_360(f.x);

% Read HYCOM data
lon0 = ncread(fhycom, 'lon');
nx0 = length(lon0);
lon0 = [lon0; lon0(1)+360];
lat0 = ncread(fhycom, 'lat');
ny0 = length(lat0);
depth0 = ncread(fhycom, 'depth');
nz0 = length(depth0);
t0 = ncread(fhycom, 'water_temp');
t0(nx0+1,:,:,:) = t0(1,:,:,:); 
s0 = ncread(fhycom, 'salinity');
s0(nx0+1,:,:,:) = s0(1,:,:,:); 
if strcmp(f.type, 'Regional')
    ix = find(lon0>=min(x360) & lon0<=max(x360));
    iy = find(lat0>=min(f.y) & lat0<=max(f.y));
    x1 = ix(1);
    nx = length(ix);
    y1 = iy(1);
    ny = length(iy);
    lon0 = lon0(ix);
    lat0 = lat0(iy);
    t0 = t0(ix,:,:,:);
    t0 = t0(:,iy,:,:);
    s0 = s0(ix,:,:,:);
    s0 = s0(:,iy,:,:);
end
[yy0, xx0] = meshgrid(lat0, lon0);
time0 = ncread(fhycom, 'time')/24 + datenum(2000,1,1);

% Interpolation
t = nan(f.node, nz0);
s = nan(f.node, nz0);
wh = interp_2d_calc_weight('BI', xx0, yy0, x360, f.y);
for iz = 1 : nz0
    disp(['Interpolating the ' num2str(iz) 'th layer of ' num2str(nz0) ' layers.'])
    % Horizontal interpolation
    t_layer = interp_2d_via_weight(t0(:,:,iz), wh);
    s_layer = interp_2d_via_weight(s0(:,:,iz), wh);
    % Fill the points that are recognize as land in HYCOM
    i_land = find(isnan(t_layer));
    i_ocean = find(~isnan(t_layer));
    k_land = knnsearch([f.x(i_ocean) f.y(i_ocean)], [f.x(i_land) f.y(i_land)], 'K', 2);
    k_land = i_ocean(k_land(:,2));
    t_layer(i_land) = t_layer(k_land);
    s_layer(i_land) = s_layer(k_land);
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

% Write initial TS output
write_initial_ts(fout, -depth0, t, s, time0);
