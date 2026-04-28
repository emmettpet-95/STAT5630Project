% Density Distributions to Orbit Errors
% Emmett Peters Spring 2026

clearvars; close all; clc;

addpath('Functions/'); addpath('Analysis')
addpath('Saved Data/'); addpath('Plots')

load('X0.mat')
load('HASDM_Grid')

%% Propagate Orbit Through Mean
clc;
tic;

time_out = HASDM_Grid.time;
time_s = seconds(time_out - time_out(1));
t0 = 0; tf = time_s(end) + 2*86400;
dt = 60*10;
tspan = t0:dt:tf;
dt_0 = time_out(1);
C_D = 2.2;

% tspan = tspan(1:1000);
% F = build_density_interpolants(HASDM_Grid);
F = 0;

[Xout, rho, r_lla] = generate_X(tspan, X0, dt_0, HASDM_Grid, time_out, C_D, F);

fprintf('Simulation Time %.1f', toc)

MeanRun.tspan = tspan;
MeanRun.tspan_dt = dt_0 + seconds(tspan);
MeanRun.X = Xout;
MeanRun.lla = r_lla;
MeanRun.lat = r_lla(1,:);
MeanRun.lon = r_lla(2,:);
MeanRun.alt = r_lla(3,:);
MeanRun.rho = rho;
% save('MeanRun.mat','MeanRun')
%% Plots
load('MeanRun.mat')
plotOrbit3D(MeanRun.X)
plotLLA(MeanRun.tspan, MeanRun.lla)
plotRho(MeanRun.tspan_dt, MeanRun.rho)

%% No Drag Run
clc;
tic;

time_out = HASDM_Grid.time;
time_s = seconds(time_out - time_out(1));
t0 = 0; tf = time_s(end)+86400*2;
dt = 60*10;
tspan = t0:dt:tf;
dt_0 = time_out(1);
scale = 0;
C_D = 2.2*scale;

% tspan = tspan(1:1000);
% F = build_density_interpolants(HASDM_Grid);
F=0;

[Xout, rho, r_lla] = generate_X(tspan, X0, dt_0, HASDM_Grid, time_out, C_D, F);

fprintf('Simulation Time %.1f', toc)

TestRun.tspan = tspan;
TestRun.tspan_dt = dt_0 + seconds(tspan);
TestRun.X = Xout;
TestRun.lla = r_lla;
TestRun.lat = r_lla(1,:);
TestRun.lon = r_lla(2,:);
TestRun.alt = r_lla(3,:);
TestRun.rho = scale*rho;
% save('NoDragRun.mat','TestRun')

%% Plots
load('NoDragRun.mat')
load('MeanRun.mat')
Xref = TestRun.X;
X = MeanRun.X;

% RIC Error
RIC_error = get_RIC(X, Xref);
Intrack = RIC_error(2,:)./1000;
figure;
plot(MeanRun.tspan/3600,Intrack)
xlabel('Time [hr]')
ylabel('In-track Error [km]')
title('In-Track Error from No Drag')
 
