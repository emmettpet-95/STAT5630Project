% Density MCMC Regression
% Emmett Peters Spring 2026

clearvars; close all; clc;

addpath('Functions/'); addpath('Analysis')
addpath('Saved Data/'); addpath('Plots')

% Get data
load('Exp2_data.mat')
%{
A = Exp2_data
A(:,1) → YYYYMMDDhhmm (timestamp)
A(:,2) → LAT (deg) -90 to 90
A(:,3) → LST (hour)
A(:,4) → RHO (kg/m^3)
A(:,5) → Dst (nT)
A(:,6) → ap (nT)
A(:,7) → f10.7
%}

alt = 400; % km

% Regression
clc;
% log(RHO) = beta0 + beta1*LAT + beta2*LST + beta3*Dst + beta4*ap + beta5*f10.7
% + epsilon(0,tau^2)

% Normalize data
A = Exp2_data;
LST_feature = cos(2*pi*(A(:,3) - 12)/24);
X = [A(:,2), LST_feature, A(:,5), A(:,6), A(:,7)];
mu_X = mean(X);
sd_X = std(X);
Xn = (X - mu_X) ./ sd_X;
y = log(A(:,4)); % log density

% Design matrix
Phi = [ones(size(Xn,1),1), Xn];   % [b0 LAT LST Dst ap F107]

% Regression
rng(12);
n_iter = 10000;
step_size = 3e-3;
n_steps = 20;
burn = 2000;
tau = std(y); % noise sd
[samples_raw, accept_rate] = run_HMC(Phi, y, tau, n_iter, step_size, n_steps);
samples = samples_raw(:,burn:end);

sample_mean = mean(samples,2);
sample_var = var(samples,0,2);
sample_ci = prctile(samples',[2.5 97.5]);

%% Plots
close all; clc
labels = {'$\beta_0$ [intercept]', '$\beta_1$ [LAT]', '$\beta_2$ [LST]', '$\beta_3$ [Dst]', '$\beta_4$ [Ap]', '$\beta_5$ [f10.7]', '$\tau$ [noise]'};
p = 6;
tspan = burn:1:size(samples_raw,2);

% Trace Plots
figure;
for i = 1:p
    subplot(p,1,i)
    plot(tspan,samples(i,:))
    ylabel(labels{i},'Interpreter','latex')
    grid on
end
xlabel('Iteration (Post Burn-in)')
sgtitle('Trace Plots')

% Posterior Dist
figure;
for i = 1:p
    subplot(2,3,i)
    histfit(samples(i,:), 30)
    hold on
    mu = mean(sample_mean(i,:));
    ci = sample_ci(:,i);
    xline(mu,'r','LineWidth',1.5)
    xline(ci(1),'m--')
    xline(ci(2),'m--')
    title(labels{i},'Interpreter','latex')
    grid on
end
sgtitle('Posterior Distributions')

% Mean Betas
figure;
bar([sample_mean(2:end)])
xticklabels(labels(2:6))
ax = gca;
ax.TickLabelInterpreter = 'latex';
ylabel('Mean \beta')
grid on
title('Posterior Means')

%% Rhat
n_chains = 3;
chains = zeros(6, n_iter, n_chains);
for c = 1:n_chains
    rng(12+c);
    chains(:,:,c) = run_HMC(Phi, y, tau, n_iter, step_size, n_steps);
end
chains = chains(:, burn:end, :);
Rhat = compute_rhat(chains);

%% Table
T = table(sample_mean, sample_var, sample_ci(1,:)', sample_ci(2,:)', Rhat, ...
    'VariableNames', {'Mean','Variance','CI_2.5','CI_97.5','Rhat'});
disp(T)