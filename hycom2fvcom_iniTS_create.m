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
%
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
t = nan(f.node, nz);
s = nan(f.node, nz);
for iz = 1 : nz
    disp(['Interpolating the ' num2str(iz) 'th layer of ' num2str(nz) ' layers.'])
    % t
    kt = ~isnan(t0(:,:,iz));
    Ft = scatteredInterpolant(xx0(kt), yy0(kt), t0(kt), 'linear', 'nearest');
    t(:,iz) = Ft(x360, f.y);
    % s
    ks = ~isnan(s0(:,:,iz));
    Fs = scatteredInterpolant(xx0(kt), yy0(kt), s0(kt), 'linear', 'nearest');
    s(:,iz) = Fs(x360, f.y);
end

% Write initial TS output
write_initial_ts(fout, -depth0, t, s, time0);
