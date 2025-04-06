%% Example: fit Pearson's data with a straight line, using covariances given by York
% Author: anna.charvatovacampbell@cmi.gov.cz
% Created: 2025-04-04
 
clear
close all
addpath("../../")

% Load data
load("PearsonYork.mat");

% Define function to be fitted as a constraint
fun = @(mu,beta) beta(1)*mu{1}-beta(2) - mu{2};


% Initial estimate of parameters and true values
beta0 = [0,0]';
mu0 = data;

% Define options
options.method = "oefpilrs2";
options.isPlot = false;
options.verbose = false;
options.maxit = 100;

% Fit the data
result = OEFPIL(data,U,fun,mu0,beta0,options);

% Check convergence
if result.iter == options.maxit
    printf("Warning: OEFPIL did not converge \n.");
end


% Print results
printf("Best estimate of parameters \n");
printf("\t Best estimate \t Uncertainty\n");
for i=1:length(beta0)
    printf("beta_%d \t %g \t %g \n", i, result.beta(i), result.ubeta(i));
end

printf("\n");

printf("True values\n");
printf("\t Best estimate \t Uncertainty\n");
for i=1:length(mu0)
    printf("mu_%d \t %g \t %g \n", i, result.mu(i), result.umu(i));
end
