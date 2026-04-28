% Density Distributions to Orbit Errors
% Emmett Peters Spring 2026

clearvars; close all; clc;

addpath('Functions/'); addpath('Analysis')
addpath('Saved Data/')

%% Get Previous Results
tic;
index_quiet = 7;
index_storm = index_quiet + 1;

[saved_data_quiet, n, p, E] = get_data(index_quiet);
[saved_data_storm, ~, ~, ~] = get_data(index_storm);

[saved_data_quiet, ~, ~, red_line_date] = get_indices(saved_data_quiet, p);
[saved_data_storm, ~, ~, ~] = get_indices(saved_data_storm, p);
red_line_date2 = saved_data_storm.tspan(1);

% Take Weighted Stats of all 100 Particles
if E == 100 

    % Weighted Average Vectors
    saved_data_quiet = get_mu2(saved_data_quiet);
    saved_data_storm = get_mu2(saved_data_storm);    
    
    % Weighted Covariance Matrices
    saved_data_quiet = get_cov2(n, p, E, saved_data_quiet);
    saved_data_storm = get_cov2(n, p, E, saved_data_storm);
    
    % Standard Deviation
    saved_data_quiet = get_sd2(n, p, saved_data_quiet);
    saved_data_storm = get_sd2(n, p, saved_data_storm);
  
    % Skewness & Kurtosis
    saved_data_quiet = get_skew_kurt2(n, p, saved_data_quiet);
    saved_data_storm = get_skew_kurt2(n, p, saved_data_storm);
    
else
    % mu_X = nan(size(UW_X_true));
    disp('NaN Weighted Means!')
    saved_data_quiet.mu_RIC = nan(size(saved_data_quiet.UW_RIC_error));
    saved_data_storm.mu_RIC = nan(size(saved_data_storm.UW_RIC_error));
    saved_data_quiet.mu_rho = nan(size(saved_data_quiet.UW_rho));
    saved_data_storm.mu_rho = nan(size(saved_data_storm.UW_rho));
end
orange = [1, 0.5, 0];

% plot_density_violin(saved_data_quiet, saved_data_storm)
% plot_I_violin(saved_data_quiet, saved_data_storm)
% plot_intrack_stats(saved_data_quiet, saved_data_storm, red_line_date, red_line_date2, orange)
% plot_density_stats(saved_data_quiet, saved_data_storm, red_line_date, red_line_date2, orange)

fprintf('Analysis Time: %.1f\n', toc)

%% Get Stats Values
mu = saved_data_storm.mu_rho;
T = length(mu);
sd = saved_data_storm.sd_rho.';
cv = sd ./ mu;
% sk = saved_data_storm.skew_rho;
% kurt = saved_data_storm.kurt_rho;

% Make Distributions
rng(1);   
N = 1000;
mu_alpha = 1;
cv1 = mean(cv(1:floor(T/2))); 
cv2 = mean(cv(floor(T/2):end))+0.1; 

% 1. Normal - small variance
sd1 = cv1 * mu_alpha;
normal_small = mu_alpha + sd1 * randn(N,1);

% 2. Normal - large variance
sd2 = cv2 * mu_alpha;
normal_large = mu_alpha + sd2 * randn(N,1);

% 3. Skewed distribution (positive skew)
cv3 = cv1;
sd3 = cv3 * mu_alpha;

alpha = 10;
skewed = skewed_CD(mu_alpha, sd3, alpha, N);

function X = skewed_CD(mu, sd, alpha, N)
    % alpha controls skew (positive = right skew)
    
    delta = alpha / sqrt(1 + alpha^2);
    
    % Generate skew-normal
    u0 = randn(N,1);
    v  = randn(N,1);
    
    X = delta .* abs(u0) + sqrt(1 - delta^2) .* v;
    
    % Normalize to mean=0, std=1
    X = (X - mean(X)) / std(X);
    
    % Scale to desired mean and std
    X = mu + sd * X;
end

% 4. Heavy-tailed distribution (high kurtosis)
cv4 = cv1;
sd4 = cv4 * mu_alpha;

% Use Student-t for heavy tails
nu = 3;  % low DOF = high kurtosis
t_samples = trnd(nu, N, 1);

% Scale to match desired std
t_samples = t_samples / std(t_samples);
kurt = mu_alpha + sd4 * t_samples;

CD_samples = normal_small;
fprintf('\nMean: %.3f\n', mean(CD_samples));
fprintf('Std: %.3f\n', std(CD_samples));
fprintf('CV: %.3f\n', std(CD_samples)/mean(CD_samples));
fprintf('Skewness: %.3f\n', skewness(CD_samples));
fprintf('Kurtosis: %.3f\n', kurtosis(CD_samples));

a_dist.normal_small = normal_small;
a_dist.normal_large = normal_large;
a_dist.skewed = skewed;
a_dist.kurt = kurt;
% save('a_dist.mat','a_dist')

plot_a_dist(a_dist)