% Density MCMC Regression
% Emmett Peters Spring 2026

clearvars; close all; clc;

addpath('Functions/'); addpath('Analysis')
addpath('Saved Data/'); addpath('Plots')

%% Read HASDM File
filename = '202405090000_202405132100_HASDM.txt';

A = readmatrix(filename, ...
    'FileType','text', ...
    'CommentStyle',{':','#'}, ...
    'Delimiter',' ', ...
    'TreatAsMissing','null');
%{
A(:,1) → YYYYMMDDhhmm (timestamp)
A(:,2) → JulianDay (NaN where "null")
A(:,3) → HTM (km) 175 to 825
A(:,4) → LAT (deg) -90 to 90
A(:,5) → LON (deg) 0 to 360
A(:,6) → LST (hour)
A(:,7) → RHO (kg/m^3)
%}

time = datetime(string(A(:,1)), 'InputFormat','yyyyMMddHH');
time_s = seconds(time - time(1));
alt  = A(:,3);
lat  = A(:,4);
lon  = A(:,5);
LST = A(:,6);
rho  = A(:,7);

[lat_unique, ~, ilat] = unique(lat);
[lon_unique, ~, ilon] = unique(lon);
[lst_unique, ~, ilst] = unique(LST);
[alt_unique, ~, ialt] = unique(alt);
[time_unique, ~, itime] = unique(time_s);

% Define limits
alt_target = 400; % km

% Logical index for rows within your ranges
idx = A(:,3) == alt_target;

% Apply index to create shortened arrays
A_short = A(idx, :);
time_short = datetime(string(A_short(:,1)), 'InputFormat','yyyyMMddHH');

HASDM_data = nan(size(A_short,1),4);
HASDM_data(:,1) = A_short(:,1); % time YYYYMMDDhhmm
HASDM_data(:,2) = A_short(:,4); % latitude
HASDM_data(:,3) = A_short(:,6); % LT
HASDM_data(:,4) = A_short(:,7); % density kg/m3

%% Read OMNI data
filename = '20240509_20240514_OMNI.txt';

A = readmatrix(filename, ...
    'FileType','text', ...
    'CommentStyle','#');
%{
A(:,1) → YYYY
A(:,2) → DOY
A(:,3) → HR
A(:,4) → Dst nT
A(:,5) → ap 
A(:,6) → f10.7
%}

year = A(:,1);
DOY  = A(:,2);
HR   = A(:,3);
mins = zeros(size(HR));   % minutes = 0

% Convert to datetime
t = datetime(year,1,1) + days(DOY-1) + hours(HR) + minutes(mins);
t_str = string(t, 'yyyyMMddHHmm');
t_num = str2double(t_str); t_num = floor(t_num / 100);

OMNI_data = nan(size(A,1),4);
OMNI_data(:,1) = t_num; % time YYYYMMDDhhmm
OMNI_data(:,2) = A(:,4); % Dst nT
OMNI_data(:,3) = A(:,5); % ap nT
OMNI_data(:,4) = A(:,6); % f10.7

%% Refine data
t_hasdm = HASDM_data(:,1);
t_omni  = OMNI_data(:,1);

% Find matching indices
[tf, loc] = ismember(t_hasdm, t_omni);

% Keep only matching HASDM rows
HASDM_filtered = HASDM_data(tf,:);

% Append corresponding OMNI data
Exp2_data = [ ...
    HASDM_filtered, ...
    OMNI_data(loc(tf), 2:4) ...
];
% save('Exp2_data.mat','Exp2_data')