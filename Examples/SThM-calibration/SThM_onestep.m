%% SThM minimalization

clear all;
pkg load statistics;
addpath("../../");

%% Set up data

caldata = dlmread("XY.dat");
x = caldata(2:end, 2);
ux = caldata(2:end, 3);
y = caldata(2:end, 4);
uy = caldata(2:end, 5);

unknown = dlmread("Y.dat");
z = unknown(2:end,2);
uz = unknown(2:end,3);

Ux = diag(ux.^2);
Uy = diag(uy.^2);
Uz = diag(uz.^2);

data = [x', y', z']';
udata = [ux', uy', uz']';

nx = length(x);
nz = length(z);

U = diag(udata.^2);


%% Set the function of parameter constraints

mu0 = data;
beta0=[[1,1,1], zeros(size(z'))]';
np =3;

fun = @(mu,beta) funSThM(mu,beta);
options.numDiffMethod = 'LynnesMoller';

%% Other options

options.tol = 1e-10;
options.method = 'oefpil';
options.criterion = 'parameterdifferences';
options.isEstimatedVariance=false;
options.maxit = 1000;

%% Run the OEFPIL algorithm

options.q = nx + nz;
result = OEFPIL(data,U,fun,mu0,beta0,options);

disp(['Number of iterations  ', num2str(result.iter)])
fprintf("\n");

disp('Parameters of calibration curve ');
fprintf("sample \t estimate \t uncertainty\n");
abc=["a","b", "c"];
for i = 1:np
    fprintf('%s %g %g\n', abc(i), result.beta(i), result.ubeta(i))
end

fprintf("\n");
disp('Conductivity of unknown samples')
fprintf("sample \t estimate \t uncertainty\n");
unknownname=["A", "B", "D", "E"];
for i = 1:nz
    fprintf('%s %g %g\n', unknownname(i), result.beta(np+i), result.ubeta(np+i))
end





