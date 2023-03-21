%==========================================================================
% HYCOM2FVCOM:
%   Create FVCOM nesting forcing from HYCOM dataset 
%
% input  :
%   fnesting  --- nesting grid file (mat)
%   ftide     --- tide constituent file (mat)
%   dir_hycom --- directory of downloaded hycom data 
%   fout      --- output nesting forcing (NetCDF)
%   fout_tide --- output nesting tide forcing (NetCDF, set '' if no need)
%   fout_hycom--- output nesting hycom forcing (NetCDF, set '' if no need)
%   time1     --- starting date and time (datenum)
%   time2     --- ending date and time (datenum)
%   dt_hycom  --- hycom file time interval (hour)
%   dt_out    --- output time interval (second)
%   lpf_T_zeta--- low-pass filter half amplitude period for elevation
%                 (hour), set nan if disabled
%   lpf_T_uv  --- low-pass filter half amplitude period for current
%                  (hour), set nan if disabled
%   zeta_adj  --- HYCOM elevation adjusted height (m)
%
% output :
%   nesting forcing file
%
% Siqi Li, SMAST
% 2022-11-30
%
% Updates:
% 2023-02-02  Siqi Li  Added two outputs for tide and hycom seperately
% 2023-02-03  Siqi Li  Removed the wrong data in HYCOM
% 2023-03-02  Siqi Li  Added low-pass filter for elevation and currents
% 2023-03-04  Siqi Li  Added an option to adjust HYCOM elevation reference
%                      level
% 2023-03-21  Siqi Li  knnsearch --> ksearch
%==========================================================================
addpath('~/tools/matFVCOM')
addpath('~/tools/t_tide')

clc
clear

%--------------------------------------------------------------------------
% Input
fnesting = '../output/gom7_nesting_grid_e.mat';
ftide = '../output/gom7_nesting_tide.mat';
dir_hycom = '../hycom1';
fout = '../output/gom7_nesting_f.nc';
fout_tide = '../output/gom7_nesting_tide_f.nc';
fout_hycom = '../output/gom7_nesting_hycom_f.nc';
fout_lpf = '../output/gom7_nesting_lpf_f.nc';
% fout = '../output/gom7_nesting_20161221_20170201.nc';
% fout_tide = '../output/gom7_nesting_tide_20161221_20170201.nc';
% fout_hycom = '../output/gom7_nesting_hycom_20161221_20170201.nc';
time1 = datenum(2016, 12, 21, 0, 0, 0);
time2 = datenum(2017, 02, 01, 0, 0, 0);
dt_hycom = 3;
dt_out = 240; 
lpf_T_zeta = 33;
lpf_T_uv = 33;
zeta_adj = 0;
%--------------------------------------------------------------------------


% Convert the dt from second to day
dt_hycom = dt_hycom / 24;
dt_out = dt_out / 3600 / 24;

% Calculate the time
t_hycom = time1 : dt_hycom : time2;
nt_hycom = length(t_hycom);
t_out = time1 : dt_out : time2;
nt_out = length(t_out);


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


% Load tide data
load(ftide);

% Initial variables
zeta3 = nan(fn.node, nt_hycom);
t3 = nan(fn.node, fn.kbm1, nt_hycom);
s3 = nan(fn.node, fn.kbm1, nt_hycom);
u3 = nan(fn.nele, fn.kbm1, nt_hycom);
v3 = nan(fn.nele, fn.kbm1, nt_hycom);

% Read the hycom data
for it = 1 : nt_hycom

    disp(['--- Interpolate HYCOM->FVCOM: ' num2str(it,'%4.4d') ' of ' num2str(nt_hycom,'%4.4d')])
    % Check the data type
    test1 = dir([dir_hycom '/hycom_' datestr(t_hycom(it), 'yyyymmdd_HHMM') '.nc']);
    test2 = dir([dir_hycom '/hycom_' datestr(t_hycom(it), 'yyyymmdd_HHMM') '_ssh.nc']);
    if ~isempty(test1)
        type = 1;
    else
        if ~isempty(test2)
            type = 2;
        else
            error('No data found.')
        end
    end

    % File names
    switch type
        case 1
            fzeta = [dir_hycom '/hycom_' datestr(t_hycom(it), 'yyyymmdd_HHMM') '.nc'];
            fts = fzeta;
            fuv = fzeta;
        case 2
            fzeta = [dir_hycom '/hycom_' datestr(t_hycom(it), 'yyyymmdd_HHMM') '_ssh.nc'];
            fts = [dir_hycom '/hycom_' datestr(t_hycom(it), 'yyyymmdd_HHMM') '_ts3z.nc'];
            fuv = [dir_hycom '/hycom_' datestr(t_hycom(it), 'yyyymmdd_HHMM') '_uv3z.nc'];
        otherwise
            error('Wrong type. Use 1 or 2.')
    end
    
    if it == 1
        % Read depth
        depth0 = ncread(fzeta, 'depth');
        nz = length(depth0);
    end

    % Read HYCOM grid
    lon0 = ncread(fzeta, 'lon');
    lon0 = calc_lon_180(lon0);
    lat0 = ncread(fzeta, 'lat');
    % Calculate interpolation weights
    wh_node = interp_2d_calc_weight('Global_BI', lon0, lat0, nesting_lon, nesting_lat);
    wh_cell = interp_2d_calc_weight('Global_BI', lon0, lat0, nesting_lonc, nesting_latc);
    wv_node = interp_vertical_calc_weight(repmat(depth0(:)',fn.node,1), fn.deplay);
    wv_cell = interp_vertical_calc_weight(repmat(depth0(:)',fn.nele,1), fn.deplayc);

    % Read the data
    zeta0 = ncread(fzeta, 'surf_el');
    t0 = ncread(fzeta, 'water_temp');
    s0 = ncread(fzeta, 'salinity');
    u0 = ncread(fzeta, 'water_u');
    v0 = ncread(fzeta, 'water_v');

    % Horizontal interpolation
    zeta1 = interp_2d_via_weight(zeta0, wh_node);
    t1 = interp_2d_via_weight(t0, wh_node);
    s1 = interp_2d_via_weight(s0, wh_node);
    u1 = interp_2d_via_weight(u0, wh_cell);
    v1 = interp_2d_via_weight(v0, wh_cell);

    % Fill the points under the bathemetry
    zeta2 = zeta1;
    t2 = fillmissing(t1, 'nearest', 2);
    s2 = fillmissing(s1, 'nearest', 2);
    u2 = fillmissing(u1, 'nearest', 2);
    v2 = fillmissing(v1, 'nearest', 2);

    % Fill the points that are recognize as land in HYCOM
    % Node
    i_land = find(isnan(t1(:,1)));
    i_ocean = find(~isnan(t1(:,1)));
    k_land = ksearch([fn.x(i_ocean) fn.y(i_ocean)], [fn.x(i_land) fn.y(i_land)], 'K', 2);
    k_land = i_ocean(k_land(:,2));
    zeta2(i_land) = zeta2(k_land);
    t2(i_land,:) = t2(k_land,:);
    s2(i_land,:) = s2(k_land,:);
    % Cell
    ic_land = find(isnan(u1(:,1)));
    ic_ocean = find(~isnan(u1(:,1)));
    kc_land = ksearch([fn.xc(ic_ocean) fn.yc(ic_ocean) ], [fn.xc(ic_land) fn.yc(ic_land)], 'K', 2);
    kc_land = ic_ocean(kc_land(:,2));
    u2(ic_land,:) = u2(kc_land,:);
    v2(ic_land,:) = v2(kc_land,:);
    
    % Vertical interpolation
    zeta3(:,it) = zeta2;
    t3(:,:,it) = interp_vertical_via_weight(t2, wv_node);
    s3(:,:,it) = interp_vertical_via_weight(s2, wv_node);
    u3(:,:,it) = interp_vertical_via_weight(u2, wv_cell);
    v3(:,:,it) = interp_vertical_via_weight(v2, wv_cell);
    
end




% Time interpolation 
% Sometimes the data could be all zeros. Such as 2017-01-03 03:00
eps = 1e-5;
t_good = squeeze(max(abs(zeta3), [], 1))>eps;
hycom_zeta = interp_time(zeta3(:,t_good), t_hycom(t_good), t_out);
t_good = squeeze(max(abs(t3-20.), [], [1 2]))>eps;
hycom_t = interp_time(t3(:,:,t_good), t_hycom(t_good), t_out);
t_good = squeeze(max(abs(s3-20.), [], [1 2]))>eps;
hycom_s = interp_time(s3(:,:,t_good), t_hycom(t_good), t_out);
t_good = squeeze(max(abs(u3), [], [1 2]))>eps;
hycom_u = interp_time(u3(:,:,t_good), t_hycom(t_good), t_out);
t_good = squeeze(max(abs(v3), [], [1 2]))>eps;
hycom_v = interp_time(v3(:,:,t_good), t_hycom(t_good), t_out);
% wt = interp_time_calc_weight(t_hycom, t_out);
% hycom_zeta = interp_time_via_weight(zeta3, wt);
% hycom_t = interp_time_via_weight(t3, wt);
% hycom_s = interp_time_via_weight(s3, wt);
% hycom_u = interp_time_via_weight(u3, wt);
% hycom_v = interp_time_via_weight(v3, wt);

% Low-pass filter
% For zeta
if isnan(lpf_T_zeta)
  hycom_zeta_lpf = hycom_zeta;
else
  disp('----Low-pass filter for elevation:')
  tmp_zeta = pl66tn(hycom_zeta, dt_out*24, lpf_T_zeta);
  hycom_zeta_lpf = tmp_zeta';
end
% For uv
if isnan(lpf_T_uv)
  hycom_u_lpf = hycom_u;
  hycom_v_lpf = hycom_v;
else
  disp('----Low-pass filter for currents:')
  hycom_u_lpf = nan(size(hycom_u));
  hycom_v_lpf = nan(size(hycom_v));
  for iz = 1 : fn.kbm1
    disp(['    Layer ' num2str(iz, '%2.2d')])
    tmp_u = squeeze(hycom_u(:,iz,:));
    tmp_v = squeeze(hycom_v(:,iz,:));
    tmp_u = pl66tn(tmp_u, dt_out*24, lpf_T_uv);
    tmp_v = pl66tn(tmp_v, dt_out*24, lpf_T_uv);
    hycom_u_lpf(:,iz,:) = tmp_u';
    hycom_v_lpf(:,iz,:) = tmp_v';
  end
end

% Adjusted HYCOM elevation
hycom_zeta_lpf = hycom_zeta_lpf + zeta_adj;

% Predict tide
% Zeta
tide_zeta = nan(fn.node, nt_out);
for i = 1 : fn.node
    if mod(i,100) == 0
        disp(['--- Predicting tide elevation: ' num2str(i, '%4.4d') ' / ' num2str(fn.node, '%4.4d')])
    end
    tide_zeta(i,:) = t_predic(t_out, tide_zeta_struct(i), ...
                              'latitude', nesting_lat(i), ...
                              'synthesis', 0);
end
% UV
tide_uv = nan(fn.nele, fn.kbm1, nt_out);
for j = 1 : fn.nele
    if mod(j,10) == 0
        disp(['--- Predicting tide velocity: ' num2str(j, '%4.4d') ' / ' num2str(fn.nele, '%4.4d')])
    end
    for iz = 1 : fn.kbm1
        tide_uv(j, iz, :) = t_predic(t_out, tide_uv_struct(j,iz), ...
                                   'latitude', nesting_latc(j), ...
                                   'synthesis', 0);
    end
end
tide_u = real(tide_uv);
tide_v = imag(tide_uv);

% Combine the tide and hycom data
out_zeta = hycom_zeta_lpf + tide_zeta;
out_u = hycom_u_lpf + tide_u;
out_v = hycom_v_lpf + tide_v;
out_t = hycom_t;
out_s = hycom_s;


% Write out the nesting file
write_nesting(fout, fn, 'Time', t_out, 'Zeta', out_zeta, ...
                                       'Temperature', out_t, ...
                                       'Salinity', out_s, ...
                                       'U', out_u, ...
                                       'V', out_v, ...
                                       'Weight_node', fn.node_weight, ...
                                       'Weight_cell', fn.cell_weight);
% Write out the nesting_tide file
if ~isempty(fout_tide)
  write_nesting(fout_tide, fn, 'Time', t_out, 'Zeta', tide_zeta, ...
                                       'U', tide_u, ...
                                       'V', tide_v);
end
% Write out the nesting_hycom file
if ~isempty(fout_hycom)
  write_nesting(fout_hycom, fn, 'Time', t_out, 'Zeta', hycom_zeta, ...
                                       'Temperature', hycom_t, ...
                                       'Salinity', hycom_s, ...
                                       'U', hycom_u, ...
                                       'V', hycom_v);
end
% Write out the low-pass filed results
if ~isempty(fout_lpf)
  write_nesting(fout_lpf, fn, 'Time', t_out, 'Zeta', hycom_zeta_lpf, ...
                                       'U', hycom_u_lpf, ...
                                       'V', hycom_v_lpf);
end


