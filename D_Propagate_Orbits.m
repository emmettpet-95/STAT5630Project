% Density Distributions to Orbit Errors
% Emmett Peters Spring 2026

clearvars; close all; clc;

addpath('Functions/'); addpath('Analysis')
addpath('Saved Data/'); addpath('Plots')

load('X0.mat')
load('HASDM_Grid')
load('CD_dist.mat')
load('a_dist.mat')
load('MeanRun.mat')
load('NoDragRun.mat')

%% Propagate Orbits Through Dist
clc;
tic;

time_out = HASDM_Grid.time;
time_s = seconds(time_out - time_out(1));
t0 = 0; tf = time_s(end)+86400*2;
dt = 60*10;
tspan = t0:dt:tf;
dt_0 = time_out(1);
N = length(tspan);
C_D = 2.2;

% CHANGE
a_array = a_dist.normal_small;
% a_array = a_dist.normal_large;
% a_array = a_dist.skewed;
% a_array = a_dist.kurt;

n_samples = 200;
a_array = a_array(1:n_samples);
A = length(a_array);

Xref = TestRun.X;
X = nan(6,N,A);
rho = nan(N,A);
r_lla = nan(3,N,A);
RIC_error = nan(3,N,A);
% F = build_density_interpolants(HASDM_Grid);
F=0;
for c = 1:A

    a = a_array(c);
    C_Dc = C_D * a;

    [Xc, unscaled_rho, r_lla(:,:,c)] = generate_X(tspan, X0, dt_0, HASDM_Grid, time_out, C_Dc, F);
    rho(:,c) = a * unscaled_rho;
    X(:,:,c) = Xc;
    RIC_error(:,:,c) = get_RIC(Xc, Xref);
    disp(c/A)

end

I_stats = compute_stats(RIC_error(2,:,:)./1000);

fprintf('Elapsed Time %.1f', toc)

Run.tspan = tspan;
Run.tspan_dt = dt_0 + seconds(tspan);
Run.X = X;
Run.lla = r_lla;
Run.rho = rho;
Run.RIC = RIC_error;
Run.Istats = I_stats;
% CHANGE
Run4 = Run;
% save('Run4.mat','Run4')

%% Plot
clc;
load('Run1.mat')
plot_intrack_stats2(Run1)
%% Plots
close all;
idx = length(Run4.tspan);
n_samples = size(Run4.X,3);
plot_a_dist(a_dist, n_samples)
plot_intrack_dist(idx)
plot_intrack_stats3()
plot_sk_k()
