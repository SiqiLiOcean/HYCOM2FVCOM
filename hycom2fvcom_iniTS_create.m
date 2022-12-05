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


% Read HYCOM data
lon0 = ncread(fhycom, 'lon');
lat0 = ncread(fhycom, 'lat');
ix = find(lon0>=min(f.x) & lon0<=max(f.x));
iy = find(lat0>=min(f.y) & lat0<=max(f.y));
x1 = ix(1);
nx = length(ix);
y1 = iy(1);
ny = length(iy);
lon0 = ncread(fhycom, 'lon', x1, nx);
lat0 = ncread(fhycom, 'lat', y1, ny);
depth0 = ncread(fhycom, 'depth');
s0 = ncread(fhycom, 'salinity', [x1 y1 1 1], [nx ny Inf 1]);
t0 = ncread(fhycom, 'water_temp', [x1 y1 1 1], [nx ny Inf 1]);
[yy0, xx0] = meshgrid(lat0, lon0);
nz = length(depth0);
time0 = ncread(fhycom, 'time')/24 + datenum(2000,1,1);


% Interpolation
t = nan(f.node, nz);
s = nan(f.node, nz);
for iz = 1 : nz
    disp(['Interpolating the ' num2str(iz) 'th layer of ' num2str(nz) ' layers.'])
    % t
    kt = ~isnan(t0(:,:,iz));
    Ft = scatteredInterpolant(xx0(kt), yy0(kt), t0(kt), 'linear', 'nearest');
    t(:,iz) = Ft(f.x, f.y);
    % s
    ks = ~isnan(s0(:,:,iz));
    Fs = scatteredInterpolant(xx0(kt), yy0(kt), s0(kt), 'linear', 'nearest');
    s(:,iz) = Fs(f.x, f.y);
end

% Write initial TS output
write_initial_ts(fout, -depth0, t, s, time0);

