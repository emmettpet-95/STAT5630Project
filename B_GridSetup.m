% Density Distributions to Orbit Errors
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
rho  = A(:,7);

%% Visualize Grid
[lat_unique, ~, ilat] = unique(lat);
[lon_unique, ~, ilon] = unique(lon);
[alt_unique, ~, ialt] = unique(alt);
[time_unique, ~, itime] = unique(time_s);

figure; hold on;
worldmap('World');
load coastlines
plotm(coastlat, coastlon, 'k')   

% Draw latitude lines
for i = 1:length(lat_unique)
    plotm(lat_unique(i)*ones(size(lon_unique)), lon_unique);
end
% Draw longitude lines
for j = 1:length(lon_unique)
    plotm(lat_unique, lon_unique(j)*ones(size(lat_unique)));
end

title(sprintf(['Lat-Lon Grid ' ...
    '\nLatitude [-90, 90] d10 (deg) \nLongitude [0, 360] d15 (deg) \nAltitude [175, 825] d25 (km) \n Time [May 9-13] d3 (hours)']));

%% Initial Simulation (No Drag)
% Earth
Req = 6378136.3e-3; %km
mu = 3.986004415e5; %km^2/s^3
J2 = 1.081874e-3;

% Orbital Elements
C_D = 2.2;
a = 400+Req; %km
ecc = 0;
inc = 45; %deg
Omega = 0;
omega = 0;
theta = 0;
oe_0 = [a, ecc, Omega, i, omega, theta];

% Initial X
X0 = oe2rv(oe_0);
% save('X0.mat','X0')

% Time
t0 = 0; tf = time_s(end);
dt = 60*1;
tspan = t0:dt:tf;
N = length(tspan);
tspan_dt = time(1) + seconds(tspan);
jds = datetime_to_jd(tspan_dt);

% Simulate
X_ND = generate_XND(tspan,X0);
plotOrbit3D(X_ND)

% Get lat and alt range
lat_ND = nan(1,N);
lon_ND = nan(1,N);
alt_ND = nan(1,N);
for k = 1:N
    r_eci = X_ND(1:3,k);
    jd = jds(k);
    [r_lla,~] = eci_to_llaecef(r_eci, jd);
    lat_ND(k) = r_lla(1);
    lon_ND(k) = r_lla(2);
    alt_ND(k) = r_lla(3);
end
r_lla_ND = [lat_ND; lon_ND; alt_ND];
plotLLA(tspan_dt, r_lla_ND)

NoDrag.tspan = tspan;
NoDrag.tspan_dt = tspan_dt;
NoDrag.X = X_ND;
NoDrag.lat = lat_ND;
NoDrag.lon = lon_ND;
NoDrag.alt = alt_ND;
% save('NoDrag.mat','NoDrag')

%% Shorten Density Array
%{
A(:,1) → YYYYMMDDhhmm (timestamp)
A(:,2) → JulianDay (NaN where "null")
A(:,3) → HTM (km) 175 to 825
A(:,4) → LAT (deg) -90 to 90
A(:,5) → LON (deg) 0 to 360
A(:,6) → LST (hour)
A(:,7) → RHO (kg/m^3)
%}

% Define limits
lat_min = -26;
lat_max = 26;
alt_min = 299;  % km
alt_max = 426;  % km
time_min = 2.024051000000000e+09;
time_max = 2.024051200000000e+09;

% Logical index for rows within your ranges
idx = (A(:,4) >= lat_min) & (A(:,4) <= lat_max) & ...
      (A(:,3) >= alt_min) & (A(:,3) <= alt_max) & ...
      (A(:,1) >= time_min) & (A(:,1) <= time_max);

% Apply index to create shortened arrays
A_short = A(idx, :);
time_short = datetime(string(A_short(:,1)), 'InputFormat','yyyyMMddHH');
lat_short = A_short(:,4);
lon_short = A_short(:,5);
alt_short = A_short(:,3);
rho_short = A_short(:,7);

%% Interpolate for finer grid
lat_step = 1;  % latitude
lon_step = 5;  % longitude
alt_step = 5;  % km altitude
time_step = 1; % mins

[Rho_out, LatGrid, LonGrid, AltGrid, time_out] = refineHasdm(A_short, lat_step, lon_step, alt_step, time_step);

HASDM_Grid.rho = Rho_out;
HASDM_Grid.lat = LatGrid;
HASDM_Grid.lon = LonGrid;
HASDM_Grid.alt = AltGrid;
HASDM_Grid.time = time_out;
% save('HASDM_Grid','HASDM_Grid')

%% Visualize finer Grid
lat_unique = -20:1:20;
lon_unique = 0:5:360;

figure; hold on;
worldmap('World');
load coastlines
plotm(coastlat, coastlon, 'k')   

% Draw latitude lines
for i = 1:length(lat_unique)
    plotm(lat_unique(i)*ones(size(lon_unique)), lon_unique);
end
% Draw longitude lines
for j = 1:length(lon_unique)
    plotm(lat_unique, lon_unique(j)*ones(size(lat_unique)));
end

title(sprintf(['Lat-Lon Grid ' ...
    '\nLatitude [-20, 20] d1 (deg) \nLongitude [0, 360] d5 (deg) \nAltitude [300, 425] d5 (km) \n Time [May 10-11] d1 (min)']));