function result = OEFPIL(data,U,fun,mu0,beta0,options)
% OEFPIL Estimates the coefficients of a nonlinear (n-variate)
% Errors-in-Variables (EIV) model specified by constraints on its parameters.
%
% We assume a measurement model: X = mu + error with fun(mu, beta) = 0,
% where X represents an N-dimensional (random) column vector of direct
% measurements, and fun represents q constraints on the model parameters,
% given by the implicit nonlinear vector function. In a more general
% context, X can be split into n subvectors. Further, beta is a
% p-dimensional column vector of model parameters that are of primary
% interest, and mu is an N-dimensional column vector representing the true
% unknown values (expectations of X), considered as the model parameters of
% secondary interest.
%
% Here, the data is represented as a (m x n)-dimensional matrix or an
% n-dimensional cell array of split (observed) subvectors x1, x2, ..., xn,
% collectively forming the N-dimensional vector of observations x = (x1',
% ..., xn')'.
%
% The covariance matrix Sigma of X = (X1',...,Xn')' is specified either by
% the uncertainty matrix U as Sigma = U, or as proportional to the
% uncertainty matrix U, Sigma = sigma2*U, where the uncertainty matrix U is
% assumed to be a known (N x N)-dimensional matrix and sigma2 is considered
% to be an unknown scalar variance parameter. U can be specified as an
% (n x n) cell array of block matrices, U{i,j} = cov(Xi, Xj) for i, j =
% 1,...,n.
%
% SYNTAX:
%    result = OEFPIL(data, U, fun, mu0, beta0, options)
%
% INPUTS
%  data    - (m x n)-dimensional matrix of measurements of the random
%            vectors X1, X2, ..., Xn. Alternatively, data is an
%            n-dimensional cell array of subvectors, i.e., data
%            = {x1, x2, ..., xn}, where xi, i = 1, ..., n are
%            (m-dimensional) column vectors.
%  U       - (N x N)-dimensional uncertainty matrix, or (n x n)-dimensional
%            cell array of block matrices U{i,j} = cov(Xi, Xj) for i, j =
%            1,...,n:
%            - If empty, the default value is U = eye(N).
%            - If U is a cell array of size (n x 1) or (1 x n), then U is a
%              block-diagonal matrix with its block matrices specified as U
%              = {U11, U22, ..., Unn}.
%            - If U is an (n x n) cell array, its block matrices are
%              specified as U = {U11, U12, ..., U1n; U21, U22, ..., U2n;
%              ... ; Un1, ...,Unn}.
%            - As U is symmetric, Uji = Uij'. If Uji is not equal to Uij',
%              we set Uji = Uij' for j < i.
%            - If any diagonal block is empty, we set Uii = eye(size(Xi)).
%            - If any off-diagonal block Uij, i>j, is empty, we set Uij =
%              zeros(size(Xi*Xj')) and Uji = zeros(size(Xj*Xi')).
%            - If any block of the cell array is a vector, it is assumed
%              that the block is a square diagonal matrix, and we set Uij =
%              diag(Uij).
%  fun     - Anonymous q-dimensional column vector function of arguments
%            mu = {mu1, mu2, ..., mun} (i.e., the elements of the (n x 1)
%            cell array) and beta (p-dimensional vector), where fun(mu,
%            beta) = 0 represents the q constraints on the model parameters
%            mu and beta.
%  mu0     - (m x n)-dimensional matrix or n-dimensional cell array of
%            m-dimensional vectors of the initial values of the parameters
%            specified in mu. Alternatively, mu0 is an n-dimensional cell
%            array of subvectors, mu0 = {mu01, ..., mu0n}. If
%            empty, we set mu0 = {x1, ..., xn}.
%  beta0   - The p-dimensional vector of the initial values of the
%            parameter vector beta0.
%  options - Structure with default values of the control parameters:
%            options.criterion = 'default';              % 'function'
%            options.criterion = 'function';             % default
%            options.criterion = 'weightedresiduals';    % alternative
%            options.criterion = 'parameterdifferences'; % alternative
%            options.maxit     = 100;
%            options.tol       = 1e-10;
%            options.delta     = eps^(1/3);
%            options.isPlot    = true;
%            options.alpha     = 0.05;
%            options.isSparse  = false;
%            options.funDiff_mu = {};      % cell array of fun derivatives
%                                          % with respect to mu{i}
%            options.funDiff_beta = [];    % fun derivatives 
%                                          % with respect to beta 
%            options.numDiffMethod = 'standard'; % standard numeric
%                                          % differentiation method
%            options.numDiffMethod = 'LynnesMoller'; % Alternative numeric
%                                          % differentiation method. Here,
%                                          % f'(x) = imag(f(x+1i*h)/h) 
%                                          % for small value of h    
%            options.h = 1e-20;            % the h value for LynnesMoller
%                                          % differentiation method
%            options.method = 'oefpil';    % default method (oefpilrs1)
%            options.method = 'oefpilrs1'; % method 1 by Radek Slesinger
%            options.method = 'oefpilrs2'; % method 2 by Radek Slesinger
%            options.method = 'oefpilvw';  % simple by Viktor Witkovsky
%            options.q = [];               % specify number of constraints
%                                          % if it is different from m. If
%                                          % empty, we set the default
%                                          % value q = m.
%            options.isEstimatedVariance = false; % If true, the algorithm
%                                          % assumes that the covariance
%                                          % matrix of measurements is
%                                          % Sigma = sigma2*U (instead of
%                                          % Sigma = U) and estimates
%                                          % the scalar parameter sigma2
%
% % EXAMPLE 1 (1-variate EIV model: Mass calibration)
%  if exist('funMassCalibration','file')
%  data  = [-9.7550e-08 -9.0350e-08  6.3200e-09 -4.7849e-07 9.4810e-08 ...
%            3.8381e-07 -3.2604e-07 -9.3040e-08  2.3331e-07 6.4100e-07 ...
%            1.2544e-04  3.0000e-07  3.6500e-02  3.6500e-02 3.0000e-02 ...
%            9.9986e-01  4.8000e-05  8.0000e-01  4.0000e-01 9.0000e-01 ...
%            1.1552e+00  5.8770e-01  8.1380e-01]';
%  ux    = [ 5.0000e-10  3.0000e-10  2.0000e-10  3.0000e-10 3.0000e-10 ...  
%            3.0000e-10  2.0000e-10  3.0000e-10  2.0000e-10 4.0000e-08 ...
%            1.0000e-09  3.0000e-08  5.0000e-04  5.0000e-04 5.0000e-04 ... 
%            5.0000e-06  2.0000e-06  1.0000e-01  1.0000e-01 1.0000e-01 ...
%            5.0000e-04  5.0000e-04  5.0000e-04]';
%  U     = diag(ux.^2);
%  fun   = @(mu,beta) funMassCalibration(mu,beta)
%  mu0   = data;
%  beta0 = [0 0 0 0]';
%  clear options
%  options.q = 9;
%  options.numDiffMethod = 'LynnesMoller';
%  options.criterion     = 'function';
%  options.method        = 'oefpilrs2';
%  options.tol           = 1e-15;
%  result = OEFPIL(data,U,fun,mu0,beta0,options);
%  else
%   error('! MISSING FILE: funMassCalibration.m')
%  end
%
% % EXAMPLE 2 (2-variate EIV model: Straight-line calibration)
%  x      = [4.0030 6.7160 9.3710 12.0530 15.2660 17.3510 ...
%           20.0360 17.3690 14.7180 12.0390 9.3760 6.6970 4.0080]';
%  y      = [0 10.1910 20.1020 30.1700 42.2300 50.0500 ...
%           60.0700 50.0800 40.1150 30.0890 20.0950 10.0700 0]';
%  m      = length(x);
%  data   = {x, y};
%  n      = length(data);
%  fun    = @(mu,beta) beta(1) + beta(2) .* mu{1} - mu{2};
%  uxA    = sqrt(0.00001444);
%  uxB    = 0.0014;
%  Ux     = uxA^2*eye(m) + uxB^2*ones(m);
%  uyA    = sqrt(0.000036);
%  uyB    = sqrt(0.000036);
%  Uy     = uyA^2*eye(m) + uyB^2*ones(m);
%  U      = {Ux, []; [] Uy};
%  mu0    = {x, y};
%  beta0  = [0;1];
%  clear options
%  options.funDiff_mu    = @(mu,beta) {beta(2).*ones(size(mu{1})), -ones(size(mu{2}))};
%  options.funDiff_beta  = @(mu,beta) [ones(size(mu{1})), mu{1}];
%  options.method        = 'oefpil';
%  options.numDiffMethod = 'LynnesMoller';
%  options.criterion     = 'parameterdifferences';
%  result = OEFPIL(data,U,fun,mu0,beta0,options);
%
% % EXAMPLE 3 (2-variate EIV model: Oliver-Phar function fit)
% % fun(mu,nu,[0.75;-0.25;1.75]))
%  x      = [ 0.2505    2.6846    5.1221    7.5628   10.0018 ...
%             0.2565    2.6858    5.1255    7.5623    9.9952 ...
%             0.2489    2.6830    5.1271    7.5603   10.0003]';
%  y      = [ 0.2398    4.9412   14.2090   27.3720   44.0513 ...
%             0.2110    4.9517   14.2306   27.3937   44.0172 ...
%             0.2303    4.9406   14.2690   27.3982   44.0611]';
%  m      = length(x);
%  data   = {x, y};
%  Ux     = 0.05^2*eye(m);
%  Uyblk  = 0.1^2*eye(m/3)+ 0.05^2*ones(m/3);
%  Uy     = blkdiag(Uyblk,Uyblk,Uyblk);
%  U      = {Ux []; [] Uy};
%  fun    = @(mu,beta) beta(1).*(mu{1}-beta(2)).^beta(3) - mu{2};
%  mu01   = x;
%  mu02   = y + 0.01*randn;
%  mu0    = {mu01 mu02};
%  beta0  = [1 0 2]';
%  clear options
%  options.numDiffMethod = 'LynnesMoller';
%  options.criterion     = 'parameterdifferences';
%  options.method        = 'oefpilrs2';
%  result = OEFPIL(data,U,fun,mu0,beta0,options);
%
% % EXAMPLE 4 (3-variate EIV model: Ellipsoid fit)
%  data = [...
%    -1.7909    1.5814   -0.4781
%     0.5736   -1.6966    0.6897
%     1.5810   -1.4075   -0.8098
%     3.5863   -0.7379    0.1183
%    -0.3867   -1.0428    1.0037
%     2.8492    1.2363    0.3700
%     3.3248   -0.6750   -0.6123
%    -1.4188    1.6040   -0.7219
%     3.6043   -0.7860   -0.2555
%     3.4583    0.9631   -0.3130
%     0.4806   -1.8653   -0.4548
%     3.5146   -0.8121    0.3495
%     0.4947   -1.4210    0.7414
%     1.7172    1.7552   -0.1043
%     0.1041    1.8793   -0.1497];
%  [m,n] = size(data);
%  N     = n*m;
%  sigma = 0.075;
%  U     = sigma^2*eye(N);
%  fun = @(mu,beta) beta(1)*mu{1}.^2 + beta(2)*mu{2}.^2 + beta(3)*mu{3}.^2 ...
%        + beta(4)*mu{1}.*mu{2} + beta(5)*mu{1}.*mu{3} ...
%        + beta(6)*mu{2}.*mu{3} + beta(7)*mu{1} + beta(8)*mu{2} ...
%        + beta(9)*mu{3} - 1;
%  mu0   = data;
%  beta0 = [1 1 1 0 0 0 0 0 0]';
%  clear options
%  options.numDiffMethod = 'LynnesMoller';
%  options.criterion     = 'parameterdifferences';
%  options.method        = 'oefpilrs2';
%  result = OEFPIL(data,U,fun,mu0,beta0,options);
%
% REFERENCES
% [1]  Charvatova Campbell, A., Slesinger, R., Klapetek, P., Chvostekova,
%      M., Hajzokova, L., Witkovsky, V. and Wimmer, G. Locally best linear
%      unbiased estimation of regression curves specified by nonlinear
%      constraints on the model parameters. AMCTMT 2023 - Advanced
%      Mathematical and Computational Tools in Metrology and Testing 2023
%      Sarajevo, Bosnia and Herzegovina, 26-28 September 2023.
% [2]  Slesinger, R., Charvatova Campbell, A., Gerslova, Z., Sindlar V.,
%      Wimmer G. (2023). OEFPIL: New method and software tool for fitting
%      nonlinear functions to correlated data with errors in variables. In
%      MEASUREMENT 2023, Smolenice, Slovakia, May 29-31, 2023, 126-129.
% [3]  Charvatova Campbell, A., Gerslova, Z., Sindlar, V., Slesinger, R.,
%      Wimmer, G. (2024). New framework for nanoindentation curve fitting
%      and measurement uncertainty estimation. Precision Engineering, 85,
%      166–173.
% [4]  Kubacek, L. (1988). Foundations of Estimation Theory. (Elsevier).
% [5]  Witkovsky, V., Wimmer, G. (2021). Polycal-MATLAB algorithm for
%      comparative polynomial calibration and its applications. In AMCTM
%      XII, 501–512.
% [6]  Koning, R., Wimmer, G., Witkovsky, V. (2014). Ellipse fitting by
%      nonlinear constraints to demodulate quadrature homodyne
%      interferometer signals and to determine the statistical uncertainty
%      of the interferometric phase. Measurement Science and Technology,
%      25(11), 115001.
% [7]  Lyness, J. N., & Moler, C. B. (1967). Numerical differentiation of
%      analytic functions. SIAM Journal on Numerical Analysis, 4(2),
%      202-210.  

% Viktor Witkovsky (witkovsky@savba.sk)
% Ver.:  18-Feb-2024 11:42:07

%% CHECK THE INPUTS AND OUTPUTS
narginchk(1, 6);
if nargin < 6, options = []; end
if nargin < 5, beta0 = []; end
if nargin < 4, mu0 = []; end
if nargin < 3, fun = []; end
if nargin < 3, U = []; end

if ~isfield(options, 'criterion')
    options.criterion = 'function';
end

if ~isfield(options, 'maxit')
    options.maxit = 100;
end

if ~isfield(options, 'tol')
    options.tol = 1e-10;
end

if ~isfield(options, 'delta')
    options.delta = eps^(1/3);
end

if ~isfield(options, 'verbose')
    options.verbose = 'true';
end

if ~isfield(options, 'isPlot')
    options.isPlot = true;
end

if ~isfield(options, 'alpha')
    options.alpha = 0.05;
end

if ~isfield(options, 'isSparse')
    options.isSparse = false;
end


if ~isfield(options, 'funDiff_mu')
    options.funDiff_mu = [];
end

if ~isfield(options, 'funDiff_nu')
    options.funDiff_nu = [];
end

if ~isfield(options, 'funDiff_beta')
    options.funDiff_beta = [];
end

if ~isfield(options, 'numDiffMethod')
    options.numDiffMethod = 'standard'; 
    % options.numDiffMethod = 'LynnesMoller'; 
    % options.numDiffMethod = 'LM'; 
end

if ~isfield(options, 'h')
    options.h = 1e-20;
end

if ~isfield(options, 'method')
    options.method = 'oefpil';
    %     options.method = 'oefpilrs1';
    %     options.method = 'oefpilrs2';
    %     options.method = 'oefpilvw';
    %     options.method = 'jacobian';
end

if ~isfield(options, 'q')
    options.q = [];
end

if ~isfield(options, 'isEstimatedVariance')
    options.isEstimatedVariance = false;
end

if iscell(data)
    n = length(data);
    m = length(data{1});
    xyz = zeros(m,n);
    for j = 1:n
        xyz(:,j) = data{j};
    end
else
    xyz = data;
    [m,n] = size(xyz);
    data = mat2cell(xyz,m,ones(1,n));
end

N   = m*n;

if isempty(options.q)
    options.q = m;
end

q = options.q;

xyzVec = xyz(:);

idm = 1:m;

if isempty(U)
    U = speye(n*m);
end

if iscell(U)
    UU = U;
    [mU,nU] = size(UU);
    if (mU == 1 && nU == 1)
        U = UU{1};
    elseif (mU == 1 && nU > 1 && nU == n)
        U = sparse(n*m,n*m);
        for j = 1:nU
            if isempty(UU{j})
                UU{j} = speye(m);
            elseif isvector(UU{j})
                UU{j} = sparse(idm,idm,UU{j});
            end
            U((j-1)*m+idm,(j-1)*m+idm) = UU{j};
        end
    elseif (mU > 1 && nU == 1 && mU == n)
        U = sparse(n*m,n*m);
        for j = 1:mU
            if isempty(UU{j})
                UU{j} = speye(m);
            elseif isvector(UU{j})
                UU{j} = sparse(idm,idm,UU{j});
            end
            U((j-1)*m+idm,(j-1)*m+idm) = UU{j};
        end
    elseif (mU > 1 && nU > 1 && mU == n && nU == n)
        U = sparse(n*m,n*m);
        for j = 1:mU
            if isempty(UU{j,j})
                UU{j,j} = speye(m);
            elseif isvector(UU{j,j})
                UU{j,j} = sparse(idm,idm,UU{j,j});
            end
            U((j-1)*m+idm,(j-1)*m+idm) = UU{j,j};
        end
        for i = 1:mU
            for j = (i+1):nU
                % The matrices in cells of symmetric UU are
                % defined by upper triangular blocks.
                % That is we require that UU{j,i} = UU{i,j}';
                UU{j,i} = UU{i,j}';
                if isempty(UU{i,j})
                    UU{i,j} = sparse(m,m);
                elseif isvector(UU{i,j})
                    UU{i,j} = sparse(idm,idm,UU{i,j});
                end
                U((i-1)*m+idm,(j-1)*m+idm) = UU{i,j};
                U((j-1)*m+idm,(i-1)*m+idm) = UU{i,j}';
            end
        end
    else
        error('OEFPIL:incorrectDimension','Error. Incorrect dimenions of U matrix.')
    end
end

if ~isempty(options.alpha)
    coverageFactor = norminv(1-options.alpha/2);
end

if options.isSparse
    U = sparse(U);
else
    U = full(U);
end

if isempty(mu0)
    mu0 = xyz;
end

if iscell(mu0)
    mu0cell = mu0;
else
    mu0cell = mat2cell(mu0,m,ones(1,n));
end

mu0 = cell2mat(mu0cell);
mu0Vec = mu0(:);
residuals = xyzVec - mu0Vec;

if isempty(beta0)
    error(['OEFPIL:incorrectDimension','Error. The starting values of the vector' ...
        ' parameter beta must be specified.'])
else
    p = length(beta0);
end

idp   = 1:p;
idB11 = [idm idm];
idB12 = [idm m+idm];
idF1  = [idm idm+m m+kron(ones(1,p),idm)];
idF2  = [idm idm kron(m+idp,ones(1,m))];
Q11 = [];
Q21 = [];
Q22 = [];
Umu = [];
umu = [];
Umubeta = [];

%% ALGORITHM
tic;
% Lower triangular matrix from Choleski decomposition of U
L = chol(U,'lower');

maxit = options.maxit;
tol   = options.tol;
crit  = 100;
iter  = 0;

%% Iterations

if any(strcmpi(options.method,{'oefpil','oefpilrs1'}))
    % OEFPILRS1 / method 1 by Radek Slesinger
    while crit > tol && iter < maxit
        iter = iter + 1;
        [B1,B2,b]  = OEFPIL_matrices(fun,mu0cell,beta0,options);
        M          = B1*U*B1';
        LM         = chol(M,'lower');
        E          = LM \ B2;
        [UE,SE,VE] = svd(E);
        F          = VE * diag(1./diag(SE));
        G          = LM' \ UE(:,1:p);
        Q21        = F*G';
        LMi        = LM \ eye(q);
        Q11        = LMi'*LMi - G*G';
        Q22        = -F*F'; %%% VW Corrected !!! Changed the sign to minus
        u          = B1*residuals + b; %%% VW Corrected !!! Changed the sing to + b
        muDelta    = residuals - U*B1'*Q11*u;
        mu0Vec     = mu0Vec + muDelta;
        mu0cell    = mat2cell(mu0Vec,m*ones(1,n),1);
        betaDelta  = -Q21*u;
        beta0      = beta0 + betaDelta;
        residuals  = xyzVec - mu0Vec;
        Lresiduals   = L\residuals;
        funcritvals  = fun(mu0cell,beta0);
        funcrit      = norm(funcritvals)/sqrt(q);
        funcritvalsL = B1*muDelta + B2*betaDelta + b;
        funcritL     = norm(funcritvalsL)/sqrt(q);
        % %%%%%%%%%%%%%%%%%%%%%
        % Space for INNER CYCLE
        % %%%%%%%%%%%%%%%%%%%%%
        if strcmpi(options.criterion,'function')
            crit  = funcrit;
        elseif strcmpi(options.criterion,'weightedresiduals')
            crit  = norm(Lresiduals)/sqrt(n*m);
        elseif strcmpi(options.criterion,'parameterdifferences')
            crit  = norm([muDelta;betaDelta]./[mu0Vec;beta0])/sqrt(n*m+p);
        else
            crit  = funcrit;
        end
    end
    Ubeta   = -Q22;
    ubeta   = sqrt(diag(Ubeta));
    Umu     = U - U*B1'*Q11*B1*U;
    umu     = sqrt(diag(Umu));
    Umb     = -U*B1'*Q21';
    Umubeta = [Umu Umb; Umb' Ubeta];
elseif strcmpi(options.method,'oefpilrs2')
    % OEFPILRS2 / method 2 by Radek Slesinger
    while crit > tol && iter < maxit
        iter = iter + 1;
        [B1,B2,b]  = OEFPIL_matrices(fun,mu0cell,beta0,options);
        M            = B1*U*B1';
        LM           = chol(M,'lower');
        E            = LM \ B2;
        [QE,RE]      = qr(E,'econ');
        REi          = RE \ eye(p);
        Q22          = -REi*REi';
        u            = B1*residuals + b; %%% VW Corrected !!! Changed the sing to + b
        v            = LM \ u;
        w            = QE'*v;
        vw           = v - QE*w;
        LMvw         = LM' \ vw;
        muDelta      = residuals - U*B1'*LMvw;
        mu0Vec       = mu0Vec + muDelta;
        mu0cell      = mat2cell(mu0Vec,m*ones(1,n),1);
        betaDelta    = -RE(1:p,:) \ w;
        beta0        = beta0 + betaDelta;
        residuals    = xyzVec - mu0Vec;
        Lresiduals   = L\residuals;
        funcritvals  = fun(mu0cell,beta0);
        funcrit      = norm(funcritvals)/sqrt(q);
        funcritvalsL = B1*muDelta + B2*betaDelta + b;
        funcritL     = norm(funcritvalsL)/sqrt(q);
        % %%%%%%%%%%%%%%%%%%%%%
        % Space for INNER CYCLE
        % %%%%%%%%%%%%%%%%%%%%%
        if strcmpi(options.criterion,'function')
            crit  = funcrit;
        elseif strcmpi(options.criterion,'weightedresiduals')
            crit  = norm(Lresiduals)/sqrt(n*m);
        elseif strcmpi(options.criterion,'parameterdifferences')
            crit  = norm([muDelta;betaDelta]./[mu0Vec;beta0])/sqrt(n*m+p);
        else
            crit  = funcrit;
        end
        if options.verbose
            Q21 = REi*QE' / LM;
            Q11 = (LM' \ (eye(q) - QE*QE')) / LM;
        end
    end
    Ubeta   = -Q22;
    ubeta   = sqrt(diag(Ubeta));
    if options.verbose
        Umu     = U - U*B1'*Q11*B1*U;
        umu     = sqrt(diag(Umu));
        Umb     = -U*B1'*Q21';
        Umubeta = [Umu Umb; Umb' Ubeta];
    end
elseif strcmpi(options.method,'oefpilvw')
    % OEFPILVW / straightforward method suggested by VW
    while crit > tol && iter < maxit
        iter = iter + 1;
        [B1,B2,b]  = OEFPIL_matrices(fun,mu0cell,beta0,options);
        z            = -(b + B1*residuals);
        B1UB1        = B1*U*B1';
        B2B1UB1z     = B2'*(B1UB1\z);
        B2B1UB1B2    = B2'*(B1UB1\B2);
        betaDelta    = B2B1UB1B2 \ B2B1UB1z;
        beta0        = beta0 + betaDelta;
        zB2betaDelta = z - B2*betaDelta;
        lambda       = B1UB1 \ zB2betaDelta;
        muDelta      = residuals + U*B1'*lambda;
        mu0Vec       = mu0Vec + muDelta;
        mu0cell      = mat2cell(mu0Vec,m*ones(1,n),1);
        residuals    = xyzVec - mu0Vec;
        Lresiduals   = L\residuals;
        funcritvals  = fun(mu0cell,beta0);
        funcrit      = norm(funcritvals)/sqrt(q);
        funcritvalsL = B1*muDelta + B2*betaDelta + b;
        funcritL     = norm(funcritvalsL)/sqrt(q);
        % %%%%%%%%%%%%%%%%%%%%%
        % Space for INNER CYCLE
        % %%%%%%%%%%%%%%%%%%%%%
        if strcmpi(options.criterion,'function')
            crit  = funcrit;
        elseif strcmpi(options.criterion,'weightedresiduals')
            crit  = norm(Lresiduals)/sqrt(n*m);
        elseif strcmpi(options.criterion,'parameterdifferences')
            crit  = norm([muDelta;betaDelta]./[mu0Vec;beta0])/sqrt(n*m+p);
        else
            crit  = funcrit;
        end
    end
    Ubeta   = B2B1UB1B2\eye(p);
    ubeta   = sqrt(diag(Ubeta));
else
    % OEFPILRS2 / method 2 by Radek Slesinger
    while crit > tol && iter < maxit
        iter = iter + 1;
        [B1,B2,b]  = OEFPIL_matrices(fun,mu0cell,beta0,options);
        M            = B1*U*B1';
        LM           = chol(M,'lower');
        E            = LM \ B2;
        [QE,RE]      = qr(E,'econ');
        REi          = RE \ eye(p);
        Q22          = -REi*REi';
        u            = B1*residuals + b; %%% VW Corrected !!! Changed the sing to + b
        v            = LM \ u;
        w            = QE'*v;
        vw           = v - QE*w;
        LMvw         = LM' \ vw;
        muDelta      = residuals - U*B1'*LMvw;
        mu0Vec       = mu0Vec + muDelta;
        mu0cell      = mat2cell(mu0Vec,m*ones(1,n),1);
        betaDelta    = -RE(1:p,:) \ w;
        beta0        = beta0 + betaDelta;
        residuals    = xyzVec - mu0Vec;
        Lresiduals   = L\residuals;
        funcritvals  = fun(mu0cell,beta0);
        funcrit      = norm(funcritvals)/sqrt(q);
        funcritvalsL = B1*muDelta + B2*betaDelta + b;
        funcritL     = norm(funcritvalsL)/sqrt(q);
        % %%%%%%%%%%%%%%%%%%%%%
        % Space for INNER CYCLE
        % %%%%%%%%%%%%%%%%%%%%%
        if strcmpi(options.criterion,'function')
            crit  = funcrit;
        elseif strcmpi(options.criterion,'weightedresiduals')
            crit  = norm(Lresiduals)/sqrt(n*m);
        elseif strcmpi(options.criterion,'parameterdifferences')
            crit  = norm([muDelta;betaDelta]./[mu0Vec;beta0])/sqrt(n*m+p);
        else
            crit  = funcrit;
        end
        if options.verbose
            Q21 = REi*QE' / LM;
            Q11 = (LM' \ (eye(q) - QE*QE')) / LM;
        end
    end
    Ubeta   = -Q22;
    ubeta   = sqrt(diag(Ubeta));
    if options.verbose
        Umu     = U - U*B1'*Q11*B1*U;
        umu     = sqrt(diag(Umu));
        Umb     = -U*B1'*Q21';
        Umubeta = [Umu Umb; Umb' Ubeta];
    end
end

% Estimated results beta and mu in different formats
beta   = beta0;
muVec  = mu0Vec;
muCell = mu0cell;
mu     = reshape(muVec,m,n);

% Estimated chiSquaredStat and the scalar variance component sigma2 
% (if we assume that Sigma = sigma^2*U) 
chiSquare =  Lresiduals'*Lresiduals;
df = q-p;
sigma2Hat = chiSquare / df;
chiSquarePval = 1-chi2cdf(chiSquare,df);

% Adjusted covariance matrix of estimators
% premultiplied by the estimated sigma2
if options.isEstimatedVariance
    Ubeta   = sigma2Hat * Ubeta;
    ubeta   = sqrt(diag(Ubeta));
    Umu     = sigma2Hat * Umu;
    umu     = sqrt(diag(Umu));
    Umubeta = sigma2Hat * Umubeta;
end

mu0 = cell2mat(mu0cell);

tictoc = toc;

if options.isPlot

    if n == 2
        figure
        plot(data{1},data{2},'*',mu0cell{1},mu0cell{2},'o')
        grid on
        xlabel('x')
        ylabel('y')
        title('EIV model: Observed vs. fitted values')
    end

    if n == 3
        figure
        plot3(data{1},data{2},data{3},'*',mu0cell{1},mu0cell{2},mu0cell{3},'o')
        grid on
        xlabel('x')
        ylabel('y')
        zlabel('z')
        title('EIV model: Observed vs. fitted values')
    end

    xx = linspace(min(xyzVec),max(xyzVec),2);
    figure
    plot(xx,xx)
    hold on
    plot(xyzVec,mu0Vec,'o')
    grid on
    xlabel('observed')
    ylabel('fitted')
    title('EIV model: Observed vs. fitted values')

    figure
    plot(residuals,'*-')
    grid on
    xlabel('index')
    ylabel('residuals')
    title('EIV model: Residuals values')
end

%% TABLES Estimated model parameters beta

TABLE_beta = table;
TABLE_beta.Properties.Description = char(fun);
TABLE_beta.ESTIMATE = beta0;
TABLE_beta.STD      = ubeta;
TABLE_beta.FACTOR   = coverageFactor*ones(size(beta0));
TABLE_beta.LOWER    = beta0 - coverageFactor*ubeta;
TABLE_beta.UPPER    = beta0 + coverageFactor*ubeta;
TABLE_beta.PVAL     = 2*normcdf(-abs(beta0./ubeta));
TABLE_beta.Properties.RowNames = string(strcat('beta_',num2str((1:p)','%-d')));

TABLE_info = table;
TABLE_info.Properties.Description = 'OEFPIL convergence';
TABLE_info.N = N;
TABLE_info.n = n;
TABLE_info.m = m;
TABLE_info.p = p;
TABLE_info.q = q;
TABLE_info.ITERATIONS  = iter;
TABLE_info.CRITERION   = crit;
TABLE_info.CONSTRAINTS_SSE = funcrit;
TABLE_info.OBJECTIVE_SSE   = Lresiduals'*Lresiduals;
%TABLE_info.FUNCCRIT_LIN = funcritL;
%TABLE_info.RSS = residuals'*residuals;
TABLE_info.ChiSquare_Stat = chiSquare;
TABLE_info.df = df;
TABLE_info.ChiSquare_Pval = chiSquarePval;
TABLE_info.sigma2Hat = sigma2Hat;

%% SHOW TABLE

if options.verbose
    disp(' ------------------------------------------------------------------------------- ')
    disp(['    OEFPIL ESTIMATION METHOD = ',char(options.method)])
    disp(['    fun = ',char(fun)])
    disp(' ------------------------------------------------------------------------------- ')
    disp(TABLE_info)
    disp(' ------------------------------------------------------------------------------- ')
    disp(TABLE_beta)
    disp(' ------------------------------------------------------------------------------- ')
end

%% Results

result.Descritpion = 'OEFPIL ESTIMATION';
result.data    = data;
result.xyz     = xyz;
result.U       = U;
result.fun     = fun;
result.beta    = beta;
result.ubeta   = ubeta;
result.Ubeta   = Ubeta;
result.mu      = mu;
result.muCell  = muCell;
result.muVec   = muVec;
result.umu     = umu;
result.Umu     = Umu;
result.Umubeta = Umubeta;
result.sigma2Hat = sigma2Hat;
result.N = N;
result.n = n;
result.m = m;
result.p = p;
result.q = q;
result.df = df;
result.chiSquare = chiSquare;
result.chiSquarePval = chiSquarePval;
result.options = options;
result.muDelta = muDelta;
result.betaDelta    = betaDelta;
result.residuals    = residuals;
result.Lresiduals   = Lresiduals;
result.funcritvals  = funcritvals;
result.funcritvalsL = funcritvalsL;
result.matrix.L  = L;
result.matrix.B1 = B1;
result.matrix.B2 = B2;
result.matrix.b  = b;
result.matrix.Q11 = Q11;
result.matrix.Q12 = Q21';
result.matrix.Q21 = Q21;
result.matrix.Q22 = Q22;
result.details.idB11 = idB11;
result.details.idB12 = idB12;
result.details.idF1  = idF1;
result.details.idF2  = idF2;
result.TABLE_beta    = TABLE_beta;
result.TABLE_INFO    = TABLE_info;
result.method        = options.method;
result.isEstimatedVariance = options.isEstimatedVariance;
result.funcritL = funcritL;
result.funcrit  = funcrit;
result.crit     = crit;
result.iter     = iter;
result.tictoc   = tictoc;

end
%% FUNCTION OEFPIL_matrices
function [B10,B20,b0] = OEFPIL_matrices(fun,mu0,beta0,options)
%OEFPIL_matrices - The required OEFPIL matrices, B1, B2, and the vector b
%  calculated from the implicit function defining the restrictions on the
%  model parameters, fun(mu,beta) = 0, computed numerically by finite
%  differences.
%
% Alternativelly, if we set 'options.isNumDiffMethodLynnesMoller = true'
% the numerical differentiation is calculated according to the method by
% Lynnes and Moller (1967). In particular, the Lynnes and Moller method
% calculates f'(x) = imag(f(x+1i*h)/h) for very small h. The default value
% is set in options.h = 1e-20.
%
% SYNTAX
%  [B10,B20,b0] = OEFPIL_matrices(fun,mu0,beta0,options)
%
% EXAMPLE 1 (Standard NumDiff method / Straight-line EIV model)
%  fun   = @(mu,beta) beta(1) + beta(2)*mu{1} - mu{2};
%  mu01  = [0.0630    0.2965    0.5321    0.7641    0.9930]';
%  mu02  = [0.0026    0.2534    0.5066    0.7558    1.0017]';
%  mu0   = {mu01, mu02};
%  beta0 = [-0.0651 1.0744]';
%  clear options
%  options.delta  = 1e-8;
%  [B1,B2,b] = OEFPIL_matrices(fun,mu0,beta0,options)
%
% EXAMPLE 2 (Lynnes and Moller NumDiff method / Straight-line EIV model)
%  fun   = @(mu,beta) beta(1) + beta(2)*mu{1} - mu{2};
%  mu01  = [0.0630    0.2965    0.5321    0.7641    0.9930]';
%  mu02  = [0.0026    0.2534    0.5066    0.7558    1.0017]';
%  mu0   = {mu01, mu02};
%  beta0 = [-0.0651 1.0744]';
%  clear options
%  options.numDiffMethod = 'LynnesMoller'
%  options.h  = 1e-50;
%  [B1,B2,b] = OEFPIL_matrices(fun,mu0,beta0,options)

% Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: '18-Feb-2024 10:33:26'

%% start gradient
narginchk(3, 4);
if nargin < 4, options = []; end

if ~isfield(options, 'delta')
    options.delta = eps^(1/3);
end

if ~isfield(options, 'isSparse')
    options.isSparse = false;
end

if ~isfield(options, 'funDiff_mu')
    options.funDiff_mu = [];
end

if ~isfield(options, 'funDiff_beta')
    options.funDiff_beta = [];
end

if ~isfield(options, 'q')
    options.q = [];
end

if ~isfield(options, 'numDiffMethod')
    options.numDiffMethod = 'standard'; 
    % options.numDiffMethod = 'LynnesMoller'; 
    % options.numDiffMethod = 'LM'; 
end

if ~isfield(options, 'h')
    options.h = 1e-20;
end

n       = length(mu0);
m       = length(mu0{1});

if isempty(options.q)
    options.q = m;
end

n_mu0   = n*m;
n_beta0 = length(beta0);
delta   = options.delta;
q       = options.q;
h       = options.h;
funDiff_mu   = options.funDiff_mu;
funDiff_beta = options.funDiff_beta;

% Function fun evaluated at given mu0 and beta0
fun_0 = fun(mu0,beta0);

% Derivatives of the function fun with respect to mu,
% evaluated at given mu0 and beta0
if isempty(funDiff_mu)
    dfun_mu0 = zeros(q,n_mu0);
    mu0v = cell2mat(mu0);
    mu0v = mu0v(:);
    % Standard method for numerical differentiation
    if any(strcmpi(options.numDiffMethod,...
            {'standard','numerical','numdiff'}))
        for i = 1:n_mu0
            mu0_minus = mu0v; mu0_minus(i) = mu0_minus(i) - delta;
            mu0_plus = mu0v; mu0_plus(i) = mu0_plus(i) + delta;
            mu0_plus = mat2cell(mu0_plus,m*ones(1,n),1);
            mu0_minus = mat2cell(mu0_minus,m*ones(1,n),1);
            dfun_mu0(:,i) = real(fun(mu0_plus,beta0) - ...
                fun(mu0_minus,beta0))/2/delta;
        end
        % Alternative (Lynne and Moller) method for numerical differentiation
    elseif any(strcmpi(options.numDiffMethod,...
            {'lynnesmoller','lm','lynnes','moller','cox','alternative'}))
        for i = 1:n_mu0
            mu0i = mu0v; mu0i(i) = mu0i(i) + 1i * h;
            mu0i = mat2cell(mu0i,m*ones(1,n),1);
            dfun_mu0(:,i) = imag(fun(mu0i,beta0)/h);
        end
        % If not specified, use the standard method
    else
        for i = 1:n_mu0
            mu0_minus = mu0v; mu0_minus(i) = mu0_minus(i) - delta;
            mu0_plus = mu0v; mu0_plus(i) = mu0_plus(i) + delta;
            mu0_plus = mat2cell(mu0_plus,m*ones(1,n),1);
            mu0_minus = mat2cell(mu0_minus,m*ones(1,n),1);
            dfun_mu0(:,i) = real(fun(mu0_plus,beta0) - ...
                fun(mu0_minus,beta0))/2/delta;
        end
    end
    % Analytically specified derivatives
else
    dfun_mu0  = funDiff_mu(mu0,beta0);
    for i = 1:n
        dfun_mu0{i} = diag(dfun_mu0{i});
    end
    dfun_mu0 = cell2mat(dfun_mu0);
end

% Derivatives of the function fun with respect to beta, evaluated at given
% mu0 and beta0
if isempty(funDiff_beta)
    dfun_beta0 = zeros(q,n_beta0);
    % Standard method for numerical differentiation
    if any(strcmpi(options.numDiffMethod,...
            {'standard','numerical','numdiff'}))
        for i = 1:n_beta0
            beta0_minus = beta0; beta0_minus(i) = beta0_minus(i) - delta;
            beta0_plus = beta0; beta0_plus(i) = beta0_plus(i) + delta;
            dfun_beta0(:,i) = real(fun(mu0,beta0_plus) - ...
                fun(mu0,beta0_minus))/2/delta;
        end
        % Alternative (Lynne and Moller) method for numerical differentiation
    elseif any(strcmpi(options.numDiffMethod,...
            {'lynnesmoller','lm','lynnes','moller','cox','alternative'}))
        for i = 1:n_beta0
            beta0i = beta0; beta0i(i) = beta0i(i) + 1i * h;
            dfun_beta0(:,i) = imag(fun(mu0,beta0i)/h);
        end
        % If not specified, use the standard method
    else
        for i = 1:n_beta0
            beta0_minus = beta0; beta0_minus(i) = beta0_minus(i) - delta;
            beta0_plus = beta0; beta0_plus(i) = beta0_plus(i) + delta;
            dfun_beta0(:,i) = real(fun(mu0,beta0_plus) - ...
                fun(mu0,beta0_minus))/2/delta;
        end
    end
    % Analytically specified derivatives
else
    dfun_beta0  = funDiff_beta(mu0,beta0);
end

% B1 = \partial fun / \partial (mu) = funDmu0
B10 = dfun_mu0;

% B2 = \partial fun / \partial (beta) = funDbeta0
B20 = dfun_beta0;

% b  = fun(mu0,nu0,beta0)
b0  = fun_0;

end