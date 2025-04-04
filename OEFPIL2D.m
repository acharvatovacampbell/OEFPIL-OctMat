function result = OEFPIL2D(x,y,U,fun,mu0,nu0,beta0,options)
%  OEFPIL2D Estimates the coefficients of a nonlinear (2-variate)
%  Errors-in-Variables (EIV) model specified by the constraints on its
%  parameters.
%
%  We assume a measurement model: X = mu + errorX, Y = nu + errorY with
%  fun(mu,nu,beta) = 0, which represents q constraints on the model
%  parameters. These constraints are given by the implicit nonlinear vector
%  function. In this context, X and Y are m-dimensional column vectors of
%  direct measurements, beta is a p-dimensional column vector of model
%  parameters that are of primary interest, and mu and nu are m-dimensional
%  column vectors representing the true unknown values (expectations of X
%  and Y), which are considered as the model parameters of secondary
%  interest.
%
%  Here x and y represent the m-dimensional vectors of measurements,
%  observed values of the random vectors X and Y.
%
%  The (N x N)-dimensional covariance matrix (uncertainty matrix), with N =
%  2*m, is assumed to be known, U = [Ux, Uxy; Uyx Uy], where Ux = cov(X,X),
%  Uxy = cov(X,Y) = Uyx', Uy = cov(Y,Y).
%
% SYNTAX:
%    result = OEFPIL2D(x,y,U,fun,mu0,nu0,beta0,options)
%
% INPUTS
%  x       - m-dimensional vector of the observed values of the input
%            quantity X = mu + errorX.
%  y       - m-dimensional vector of the observed values of the input
%            quantity Y = mu + errorY.
%  U       - (N x N)-dimensional uncertainty matrix, N = 2*m, where
%            U = [Ux, Uxy; Uyx, Uy], or is defined as a cell structure U =
%            {Ux, Uy, Uxy}, where, Ux - the (m x m)-dimensional uncertainty
%            matrix of the quantity X, Uy - the (m x m)-dimensional
%            uncertainty matrix of the quantity Y, Uxy- the (m x m)
%            dimensional covariance matrix of the quantities X and Y. If
%            empty, default value is U = {eye(m), eye(m), zeros{m}}, where
%            m = length(x).
%  fun     - the anonymous q-dimensional column vector function of
%            arguments mu (m-dimensional column vector), nu (m-dimensional
%            column vector), and beta (p-dimensional vector),  where
%            fun(mu,nu,beta) = 0 represent the q constraints on the model
%            parameters. If empty, default value is a (straight-line)
%            function, fun = @(mu,nu,beta) beta(1) + beta(2)*mu - nu.
%  mu0     - the m-dimensional vector of the initial values of the
%            parameter vector mu.
%  nu0     - the m-dimensional vector of the initial values of the
%            parameter vector nu.
%  beta0   - the p-dimensional vector of the initial values of the
%            parameter vector nbeta.
%  options - structure with the following default values of the control
%            parameters:
%            options.criterion = 'default';              % 'function'
%            options.criterion = 'function';             % default
%            options.criterion = 'weightedresiduals';    % alternative
%            options.criterion = 'parameterdifferences'; % alternative
%            options.maxit = 100;
%            options.tol = 1e-10;
%            options.delta = eps^(1/3);
%            options.isPlot = true;
%            options.alpha = 0.05;
%            options.isSparse = false;
%            options.funDiff_mu = [];      % fun derivatives 
%                                          % with respect to mu
%            options.funDiff_nu = [];      % fun derivatives 
%                                          % with respect to nu
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
%            options.method = 'jacobian';  % method based on linearized
%                                          % model with Jacobian
%            options.q = [];               % specify nymber of constraints
%                                          % if it is different from m. If
%                                          % empty we set the default value
%                                          % q = m.
%            options.isEstimatedVariance = false; % If true, the algorithm
%                                          % assumes that the covariance
%                                          % matrix of measurements is
%                                          % Sigma = sigma2*U (instead of
%                                          % Sigma = U) and estimates
%                                          % the scalar parameter sigma2
%
% EXAMPLE 1 (Straight-line calibration)
%  x      = [4.0030 6.7160 9.3710 12.0530 15.2660 17.3510 ...
%           20.0360 17.3690 14.7180 12.0390 9.3760 6.6970 4.0080]';
%  y      = [0 10.1910 20.1020 30.1700 42.2300 50.0500 ...
%           60.0700 50.0800 40.1150 30.0890 20.0950 10.0700 0]';
%  fun    = @(mu,nu,beta) beta(1) + beta(2) .* mu - nu;
%  m      = length(x);
%  uxA    = sqrt(0.00001444);
%  uxB    = 0.0014;
%  Ux     = uxA^2*eye(m) + uxB^2*ones(m);
%  uyA    = sqrt(0.000036);
%  uyB    = sqrt(0.000036);
%  Uy     = uyA^2*eye(m) + uyB^2*ones(m);
%  Uxy    = [];
%  mu0    = x;
%  nu0    = y;
%  beta0  = [0;1];
%  clear options
%  options.funDiff_mu   = @(mu,nu,beta) beta(2).*ones(size(mu));
%  options.funDiff_nu   = @(mu,nu,beta) -ones(size(nu));
%  options.funDiff_beta = @(mu,nu,beta) [ones(size(mu)), mu];
%  options.method       = 'oefpil';
%  options.criterion    = 'parameterdifferences';
%  result = OEFPIL2D(x,y,{Ux,Uy,Uxy},fun,mu0,nu0,beta0,options);
%
% EXAMPLE 2 (Oliver-Phar function fit / fun(mu,nu,[0.75;-0.25;1.75]))
%  x     = [ 0.2505    2.6846    5.1221    7.5628   10.0018 ...
%            0.2565    2.6858    5.1255    7.5623    9.9952 ...
%            0.2489    2.6830    5.1271    7.5603   10.0003]';
%  y     = [ 0.2398    4.9412   14.2090   27.3720   44.0513 ...
%            0.2110    4.9517   14.2306   27.3937   44.0172 ...
%            0.2303    4.9406   14.2690   27.3982   44.0611]';
%  k     = 5;
%  m     = 15;
%  Ux    = 0.05^2*eye(m);
%  Uyblk = 0.1^2*eye(k)+ 0.05^2*ones(k);
%  Uy    = blkdiag(Uyblk,Uyblk,Uyblk);
%  Uxy   = zeros(m);
%  U     = [Ux Uxy; Uxy' Uy];
%  fun   = @(mu,nu,beta) beta(1).*(mu-beta(2)).^beta(3) - nu;
%  mu0   = x;
%  nu0   = y + 0.01*randn;
%  beta0 = [1 0 2]';
%  clear options
%  options.method       = 'jacobian';
%  options.funDiff_mu   = @(mu,nu,beta) beta(1)*beta(3) ...
%                          *(mu - beta(2)).^(beta(3)-1);
%  options.funDiff_nu   = @(mu,nu,beta) - ones(size(nu));
%  options.funDiff_beta = @(mu,nu,beta) [(mu-beta(2)).^beta(3), ...
%                       -beta(1)*beta(3)*(mu-beta(2)).^(beta(3)-1), ...
%                        beta(1).*(mu-beta(2)).^beta(3).*log(mu-beta(2))];
%  result = OEFPIL2D(x,y,U,fun,mu0,nu0,beta0,options);
%
% EXAMPLE 3 (Flow Meter Calibration / Fractional Polynomial Calibration)
%  x   = [3.1100e+03 2.5000e+03 1.5000e+03 7.4000e+02 4.2200e+02 ...
%         1.2000e+02 2.5000e+01 1.0000e+01 7.0000e+00 5.5000e+00 ...
%         3.2500e+00]';
%  y   = [-1.4040e+00 -1.1500e+00 -6.0900e-01 2.6200e-01 9.8000e-01 ...
%         1.9290e+00 1.3400e+00 -2.5400e-01 -1.3800e+00 -2.4220e+00 ...
%         -5.2990e+00]';
%  m   = 11;
%  Ux1 = [9.7000e+00 3.1000e+00 1.1000e+00 2.7000e-01 8.7000e-02 7.1000e-03;
%         3.1000e+00 6.3000e+00 1.1000e+00 2.7000e-01 8.7000e-02 7.1000e-03;
%         1.1000e+00 1.1000e+00 2.3000e+00 2.7000e-01 8.7000e-02 7.1000e-03;
%         2.7000e-01 2.7000e-01 2.7000e-01 5.5000e-01 8.7000e-02 7.1000e-03;
%         8.7000e-02 8.7000e-02 8.7000e-02 8.7000e-02 1.8000e-01 7.1000e-03;
%         7.1000e-03 7.1000e-03 7.1000e-03 7.1000e-03 7.1000e-03 1.4000e-02];
%  Ux2 = [2.5000e-03 2.0000e-04 9.6000e-05 5.9000e-05 2.1000e-05;
%         2.0000e-04 4.0000e-04 9.6000e-05 5.9000e-05 2.1000e-05;
%         9.6000e-05 9.6000e-05 2.0000e-04 5.9000e-05 2.1000e-05;
%         5.9000e-05 5.9000e-05 5.9000e-05 1.2000e-04 2.1000e-05;
%         2.1000e-05 2.1000e-05 2.1000e-05 2.1000e-05 4.2000e-05];
%  Ux  = blkdiag(Ux1,Ux2);
%  Uy  = diag([1.0240e-03 3.7210e-03 2.6010e-03 1.7640e-03 3.7210e-03 ...
%              4.0960e-03 3.3640e-03 3.2490e-03 6.0840e-03 6.2410e-03 ...
%              1.0609e-02]);
%  Uxy = zeros(m);
%  fun = @(mu,nu,beta) beta(1)*mu.^(-1) + beta(2)*mu.^(-0.5) + ...
%        beta(3) + beta(4)*mu.^(0.5) + beta(5)*mu - nu;
%  mu0 = x;
%  nu0 = y;
%  beta0 = [0; 0; 0; 0; 1];
%  clear options
%  options.numDiffMethod = 'LynnesMoller'; 
%  options.criterion     = 'parameterdifferences';
%  options.method        = 'oefpilrs2';
%  options.tol           = 1e-12;
%  result = OEFPIL2D(x,y,{Ux,Uy,Uxy},fun,mu0,nu0,beta0,options);
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
narginchk(2, 8);
if nargin < 8, options = []; end
if nargin < 7, beta0 = []; end
if nargin < 6, nu0 = []; end
if nargin < 5, mu0 = []; end
if nargin < 4, fun = []; end
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

x   = x(:);
y   = y(:);
n   = 2;
m   = length(x);
N   = m*n;
idm = 1:m;

if isempty(options.q)
    options.q = m;
end
q = options.q;

if isempty(U)
    U = speye(N);
else
    if iscell(U)
        if isempty(U{1})
            U{1} = speye(m);
        else
            if isvector(U{1})
                U{1} = sparse(idm,idm,U{1});
            end
        end
        if isempty(U{2})
            U{2} = speye(m);
        else
            if isvector(U{2})
                U{2} = sparse(idm,idm,U{2});
            end
        end
        if isempty(U{3})
            U{3} = sparse(m,m);
        else
            if isvector(U{3})
                U{3} = sparse(idm,idm,U{3});
            end
        end
        U = [U{1} U{3}; U{3} U{2}];
    end
end

if ~isempty(options.alpha)
    coverageFactor = norminv(1-options.alpha/2);
end

if options.isSparse
    U = sparse(U);
end

if isempty(fun)
    fun = @(mu,nu,beta) beta(1) + beta(2) .* mu - nu;
    X = [ones(m,1) x];
    beta0 = (X'*X) \ (X'*y);
end

if isempty(mu0)
    mu0 = x;
end

if isempty(nu0)
    nu0 = y;
end

if isempty(beta0)
    error(['Error. The starting values of the vector' ...
        ' parameter beta must be specified'])
else
    p = length(beta0);
end

idp   = 1:p;
idB11 = [idm idm];
idB12 = [idm m+idm];
idF1  = [idm idm+m m+kron(ones(1,p),idm)];
idF2  = [idm idm kron(m+idp,ones(1,m))];
idmp  = [idm,2*m+idp];
Q11 = [];
Q21 = [];
Q22 = [];
Umunu = [];
umunu = [];
Umubeta = [];
Umunubeta = [];

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
        [B1,B2,b,J] = OEFPIL_matrices(fun,mu0,nu0,beta0,options);
        xyDelta = [x - mu0;y - nu0];
        M       = B1*U*B1';
        LM      = chol(M,'lower');
        E       = LM \ B2;
        [UE,SE,VE] = svd(E);
        F   = VE * diag(1./diag(SE));
        G   = LM' \ UE(:,1:p);
        Q21 = F*G';
        LMi = LM \ eye(m);
        Q11 = LMi'*LMi - G*G';
        Q22 = -F*F'; %%% VW Corrected !!! Changed the sign to minus
        u   = B1*xyDelta + b; %%% VW Corrected !!! Changed the sing to + b
        munuDelta = xyDelta - U*B1'*Q11*u;
        muDelta = munuDelta(idm);
        nuDelta = munuDelta(m+idm);
        betaDelta = -Q21*u;
        mu0   = mu0 + muDelta;
        nu0   = nu0 + nuDelta;
        beta0 = beta0 + betaDelta;
        xResiduals = x - mu0;
        yResiduals = y - nu0;
        Lresiduals  = L\[xResiduals;yResiduals];
        funcritvals  = fun(mu0,nu0,beta0);
        funcrit      = norm(funcritvals)/sqrt(q);
        funcritvalsL = B1*munuDelta + B2*betaDelta + b;
        funcritL     = norm(funcritvalsL)/sqrt(q);
        if strcmpi(options.criterion,'function')
            crit  = funcrit;
        elseif strcmpi(options.criterion,'weightedresiduals')
            crit  = norm(Lresiduals)/sqrt(N);
        elseif strcmpi(options.criterion,'parameterdifferences')
            crit  = norm([muDelta;nuDelta;betaDelta]./[mu0;nu0;beta0])/sqrt(N+p);
        else
            crit  = funcrit;
        end
    end
    Ubeta   = -Q22;
    ubeta   = sqrt(diag(Ubeta));
    Umunu   = U - U*B1'*Q11*B1*U;
    umunu   = sqrt(diag(Umunu));
    Umb     = -U*B1'*Q21';
    Umunubeta = [Umunu Umb; Umb' Ubeta];
    Umubeta   = Umunubeta(idmp,idmp);
    %elseif strcmpi(options.method,'oefpil2')
elseif strcmpi(options.method,'oefpilrs2')
    % OEFPILRS2 / method 2 by Radek Slesinger
    while crit > tol && iter < maxit
        iter = iter + 1;
        [B1,B2,b,J] = OEFPIL_matrices(fun,mu0,nu0,beta0,options);
        xyDelta  = [x - mu0;y - nu0];
        M            = B1*U*B1';
        LM           = chol(M,'lower');
        E            = LM \ B2;
        [QE,RE]      = qr(E,'econ');
        REi          = RE \ eye(p);
        Q22          = -REi*REi';
        u            = B1*xyDelta + b; %%% VW Corrected !!! Changed the sing to + b
        v            = LM \ u;
        w            = QE'*v;
        vw           = v - QE*w;
        LMvw         = LM' \ vw;
        munuDelta    = xyDelta - U*B1'*LMvw;
        muDelta      = munuDelta(idm);
        nuDelta      = munuDelta(m+idm);
        mu0          = mu0 + muDelta;
        nu0          = nu0 + nuDelta;
        betaDelta    = -RE(1:p,:) \ w;
        beta0        = beta0 + betaDelta;
        xResiduals   = x - mu0;
        yResiduals   = y - nu0;
        Lresiduals = L\[xResiduals;yResiduals];
        funcritvals  = fun(mu0,nu0,beta0);
        funcrit      = norm(funcritvals)/sqrt(q);
        funcritvalsL = B1*munuDelta + B2*betaDelta + b;
        funcritL     = norm(funcritvalsL)/sqrt(q);
        if strcmpi(options.criterion,'function')
            crit  = funcrit;
        elseif strcmpi(options.criterion,'weightedresiduals')
            crit  = norm(Lresiduals)/sqrt(N);
        elseif strcmpi(options.criterion,'parameterdifferences')
            crit  = norm([muDelta;nuDelta;betaDelta]./[mu0;nu0;beta0])/sqrt(N+p);
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
        Umunu   = U - U*B1'*Q11*B1*U;
        umunu   = sqrt(diag(Umunu));
        Umb     = -U*B1'*Q21';
        Umunubeta = [Umunu Umb; Umb' Ubeta];
        Umubeta   = Umunubeta(idmp,idmp);
    end
elseif strcmpi(options.method,'jacobian')
    % method based on locally linearized model by using Jacobian matrix
    % the method is useful if there is explicite functional relation
    % between nu amd mu: nu = g(mu,beta), hence the implicite condition
    % f(mu,nu,beta) = 0 is given by f(mu,nu,beta) = g(mu,beta) - nu = 0.
    % JACOBIAN / locally linearized model by VW
    while crit > tol && iter < maxit
        iter = iter + 1;
        [B1,B2,b,J] = OEFPIL_matrices(fun,mu0,nu0,beta0,options);
        F    = -J;
        FL   = L\F;
        FLLF = FL' * FL;
        LXY  = L\[x - mu0;y - nu0];
        mubetaDelta = FLLF \ (FL' * LXY);
        muDelta = mubetaDelta(idm);
        mu0   = mu0 + muDelta;
        betaDelta = mubetaDelta(m+idp);
        beta0 = beta0 + betaDelta;
        funcritvals  = fun(mu0,nu0,beta0);
        nuDelta = funcritvals;
        nu0   = nu0 + nuDelta;
        funcrit      = norm(funcritvals)/sqrt(q);
        funcritvalsL = B1*[muDelta;nuDelta] + B2*betaDelta + b;
        funcritL     = norm(funcritvalsL)/sqrt(q);
        xResiduals = x - mu0;
        yResiduals = y - nu0;
        Lresiduals  = L\[xResiduals;yResiduals];
        if strcmpi(options.criterion,'function')
            crit  = funcrit;
        elseif strcmpi(options.criterion,'weightedresiduals')
            crit  = max(abs(mubetaDelta./[mu0;beta0]));
        elseif strcmpi(options.criterion,'parameterdifferences')
            crit  = norm(mubetaDelta./[mu0;beta0]);
        else
            crit  = funcrit;
        end
    end
    Umubeta = FLLF \ eye(m+p);
    Ubeta   = Umubeta(m+idp,m+idp);
    ubeta   = sqrt(diag(Ubeta));
elseif strcmpi(options.method,'oefpilvw')
    % OEFPILVW / straightforward method suggested by VW
    while crit > tol && iter < maxit
        iter = iter + 1;
        [B1,B2,b,J]  = OEFPIL_matrices(fun,mu0,nu0,beta0,options);
        xyDelta      = [x - mu0;y - nu0];
        z            = -(b + B1*xyDelta);
        B1UB1        = B1*U*B1';
        B2B1UB1z     = B2'*(B1UB1\z);
        B2B1UB1B2    = B2'*(B1UB1\B2);
        betaDelta    = B2B1UB1B2 \ B2B1UB1z;
        beta0        = beta0 + betaDelta;
        zB2betaDelta = z - B2*betaDelta;
        lambda       = B1UB1 \ zB2betaDelta;
        munuDelta    = xyDelta + U*B1'*lambda;
        muDelta      = munuDelta(idm);
        nuDelta      = munuDelta(m+idm);
        mu0          = mu0 + muDelta;
        nu0          = nu0 + nuDelta;
        xResiduals   = x - mu0;
        yResiduals   = y - nu0;
        Lresiduals = L\[xResiduals;yResiduals];
        funcritvals  = fun(mu0,nu0,beta0);
        funcrit      = norm(funcritvals)/sqrt(q);
        funcritvalsL = B1*munuDelta + B2*betaDelta + b;
        funcritL     = norm(funcritvalsL)/sqrt(q);
        if strcmpi(options.criterion,'function')
            crit  = funcrit;
        elseif strcmpi(options.criterion,'weightedresiduals')
            crit  = norm(Lresiduals)/sqrt(N);
        elseif strcmpi(options.criterion,'parameterdifferences')
            crit  = norm([muDelta;nuDelta;betaDelta]./[mu0;nu0;beta0])/sqrt(N+p);
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
        [B1,B2,b,J] = OEFPIL_matrices(fun,mu0,nu0,beta0,options);
        xyDelta  = [x - mu0;y - nu0];
        M            = B1*U*B1';
        LM           = chol(M,'lower');
        E            = LM \ B2;
        [QE,RE]      = qr(E,'econ');
        REi          = RE \ eye(p);
        Q22          = -REi*REi';
        u            = B1*xyDelta + b; %%% VW Corrected !!! Changed the sing to + b
        v            = LM \ u;
        w            = QE'*v;
        vw           = v - QE*w;
        LMvw         = LM' \ vw;
        munuDelta    = xyDelta - U*B1'*LMvw;
        muDelta      = munuDelta(idm);
        nuDelta      = munuDelta(m+idm);
        mu0          = mu0 + muDelta;
        nu0          = nu0 + nuDelta;
        betaDelta    = -RE(1:p,:) \ w;
        beta0        = beta0 + betaDelta;
        xResiduals   = x - mu0;
        yResiduals   = y - nu0;
        Lresiduals   = L\[xResiduals;yResiduals];
        funcritvals  = fun(mu0,nu0,beta0);
        funcrit      = norm(funcritvals);
        funcritvalsL = B1*munuDelta + B2*betaDelta + b;
        funcritL     = norm(funcritvalsL);
        if strcmpi(options.criterion,'function')
            crit  = funcrit;
        elseif strcmpi(options.criterion,'weightedresiduals')
            crit  = norm(Lresiduals)/sqrt(2*m);
        elseif strcmpi(options.criterion,'parameterdifferences')
            crit  = norm([muDelta;nuDelta;betaDelta]./[mu0;nu0;beta0])/sqrt(N+p);
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
        Umunu   = U - U*B1'*Q11*B1*U;
        umunu   = sqrt(diag(Umunu));
        Umb     = -U*B1'*Q21';
        Umunubeta = [Umunu Umb; Umb' Ubeta];
        Umubeta   = Umunubeta(idmp,idmp);
    end
end

% Estimated results beta and mu in different formats
beta    = beta0;
mu      = mu0;
nu      = nu0;
munu    = [mu nu];
munuVec = munu(:);
munuCell = {mu nu};

% Estimated chiSquaredStat and the scalar variance component sigma2 
% (if we assume that Sigma = sigma^2*U) 
residuals = [xResiduals;yResiduals];
chiSquare =  Lresiduals'*Lresiduals;
df = q-p;
sigma2Hat = chiSquare / df;
chiSquarePval = 1-chi2cdf(chiSquare,df);

% Adjusted covariance matrix of estimators
% premultiplied by the estimated sigma2
if options.isEstimatedVariance
    Ubeta   = sigma2Hat * Ubeta;
    ubeta   = sqrt(diag(Ubeta));
    Umunu   = sigma2Hat * Umunu;
    umunu   = sqrt(diag(Umunu));
    Umunubeta = sigma2Hat * Umunubeta;
end

tictoc = toc;

if options.isPlot
    figure
    plot(x,y,'*',mu0,nu0,'o')
    grid on
    xlabel('x')
    ylabel('y')
    title('EIV model: Observed vs. fitted values')
    %
    %     figure
    %     plot([xResiduals;yResiduals],'*-')
    %     grid on
    %     xlabel('index')
    %     ylabel('residuals')
    %     title('EIV model: Residuals values')
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
%TABLE_info.xRSS = xResiduals'*xResiduals;
%TABLE_info.yRSS = yResiduals'*yResiduals;
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
result.x    = x;
result.y    = y;
result.xy   = [x y];
result.U    = U;
result.fun  = fun;
result.mu   = mu;
result.nu   = nu;
result.munu = munu;
result.munuCell  = munuCell;
result.munuVec   = munuVec;
result.beta      = beta;
result.ubeta     = ubeta;
result.Ubeta     = Ubeta;
result.umunu     = umunu;
result.Umunu     = Umunu;
result.Umunubeta = Umunubeta;
result.Umubeta   = Umubeta;
result.sigma2    = sigma2Hat;
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
result.nuDelta = nuDelta;
result.betaDelta    = betaDelta;
result.xResiduals   = xResiduals;
result.yResiduals   = yResiduals;
result.residuals    = residuals;
result.Lresiduals   = Lresiduals;
result.funcritvals  = funcritvals;
result.funcritvalsL = funcritvalsL;
result.matrix.L  = L;
result.matrix.B1 = B1;
result.matrix.B2 = B2;
result.matrix.b  = b;
result.matrix.J  = J;
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
function [B10,B20,b0,J0] = OEFPIL_matrices(fun,mu0,nu0,beta0,options)
%OEFPIL_matrices - The required OEFPIL matrices, B1, B2, and the vector b
%  calculated from the implicit function defining the restrictions on the
%  model parameters, fun(mu,nu,beta) = 0, computed numerically by finite
%  differences.
%
%  Alternativelly, if we set 'options.isNumDiffMethodLynnesMoller = true'
%  the numerical differentiation is calculated according to the method by
%  Lynnes and Moller (1967). In particular, the Lynnes and Moller method
%  calculates f'(x) = imag(f(x+1i*h)/h) for very small h. The default value
%  is set in options.h = 1e-20.
%
%  Moreover, if required the Jacobian matrix J of the residual vector [x -
%  mu; y - nu ] is evaluated with respect to the parameters mu and beta. J
%  is valid Jacobian matrix if the restriction can be expressed in explicit
%  form, as nu = g(mu,beta). So, here we assume that the implicit form of
%  the function fun is of the following form, fun(mu,nu,beta) = g(mu,beta)
%  -nu.
%
% SYNTAX
%  [B10,B20,b0,J0] = OEFPIL_matrices(fun,mu0,nu0,beta0,options)
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
narginchk(4, 5);
if nargin < 5, options = []; end

if ~isfield(options, 'delta')
    options.delta = eps^(1/3);
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

n_mu0   = length(mu0);
n_nu0   = length(nu0);
n_beta0 = length(beta0);

delta = options.delta;
h       = options.h;
funDiff_mu = options.funDiff_mu;
funDiff_nu = options.funDiff_nu;
funDiff_beta = options.funDiff_beta;

% Function fun evaluated at given mu0, nu0 and beta0
fun_0 = fun(mu0,nu0,beta0);

% Derivatives of the function fun with respect to mu, evaluated at given
% mu0, nu0 and beta0
if isempty(funDiff_mu)
    dfun_mu0 = zeros(n_mu0,n_mu0);
    % Standard method for numerical differentiation
    if any(strcmpi(options.numDiffMethod,...
            {'standard','numerical','numdiff'}))
        for i = 1:n_mu0
            mu0_minus = mu0; mu0_minus(i) = mu0_minus(i) - delta;
            mu0_plus = mu0; mu0_plus(i) = mu0_plus(i) + delta;
            dfun_mu0(:,i) = real(fun(mu0_plus,nu0,beta0) - ...
                fun(mu0_minus,nu0,beta0))/2/delta;
        end
        % Alternative (Lynne and Moller) method for numerical differentiation
    elseif any(strcmpi(options.numDiffMethod,...
            {'lynnesmoller','lm','lynnes','moller','cox','alternative'}))
        for i = 1:n_mu0
            mu0i = mu0; mu0i(i) = mu0i(i) + 1i * h;
            dfun_mu0(:,i) = imag(fun(mu0i,nu0,beta0)/h);
        end
        % If not specified, use the standard method
    else
        for i = 1:n_mu0
            mu0_minus = mu0; mu0_minus(i) = mu0_minus(i) - delta;
            mu0_plus = mu0; mu0_plus(i) = mu0_plus(i) + delta;
            dfun_mu0(:,i) = real(fun(mu0_plus,nu0,beta0) - ...
                fun(mu0_minus,nu0,beta0))/2/delta;
        end
    end
else
    dfun_mu0 = diag(funDiff_mu(mu0,nu0,beta0));
end

% Derivatives of the function fun with respect to nu, evaluated at given
% mu0, nu0 and beta0
if isempty(funDiff_nu)
    dfun_nu0 = zeros(n_mu0,n_nu0);
    % Standard method for numerical differentiation
    if any(strcmpi(options.numDiffMethod,...
            {'standard','numerical','numdiff'}))
        for i = 1:n_nu0
            nu0_minus = nu0; nu0_minus(i) = nu0_minus(i) - delta;
            nu0_plus = nu0; nu0_plus(i) = nu0_plus(i) + delta;
            dfun_nu0(:,i) = real(fun(mu0,nu0_plus,beta0) - ...
                fun(mu0,nu0_minus,beta0))/2/delta;
        end
    elseif any(strcmpi(options.numDiffMethod,...
            {'lynnesmoller','lm','lynnes','moller','cox','alternative'}))
        for i = 1:n_nu0
            nu0i = nu0; nu0i(i) = nu0i(i) + 1i * h;
            dfun_nu0(:,i) = imag(fun(mu0,nu0i,beta0)/h);
        end
        % If not specified, use the standard method
    else
        for i = 1:n_nu0
            nu0_minus = nu0; nu0_minus(i) = nu0_minus(i) - delta;
            nu0_plus = nu0; nu0_plus(i) = nu0_plus(i) + delta;
            dfun_nu0(:,i) = real(fun(mu0,nu0_plus,beta0) - ...
                fun(mu0,nu0_minus,beta0))/2/delta;
        end
    end
else
    dfun_nu0 = diag(funDiff_nu(mu0,nu0,beta0));
end

% Derivatives of the function fun with respect to beta, evaluated at given
% mu0, nu0 and beta0
if isempty(funDiff_beta)
    dfun_beta0 = zeros(n_mu0,n_beta0);
    % Standard method for numerical differentiation
    if any(strcmpi(options.numDiffMethod,...
            {'standard','numerical','numdiff'}))
        for i = 1:n_beta0
            beta0_minus = beta0; beta0_minus(i) = beta0_minus(i) - delta;
            beta0_plus = beta0; beta0_plus(i) = beta0_plus(i) + delta;
            dfun_beta0(:,i) = real(fun(mu0,nu0,beta0_plus) - ...
                fun(mu0,nu0,beta0_minus))/2/delta;
        end
    elseif any(strcmpi(options.numDiffMethod,...
            {'lynnesmoller','lm','lynnes','moller','cox','alternative'}))
        for i = 1:n_beta0
            beta0i = beta0; beta0i(i) = beta0i(i) + 1i * h;
            dfun_beta0(:,i) = imag(fun(mu0,nu0,beta0i)/h);
        end
        % If not specified, use the standard method
    else
        for i = 1:n_beta0
            beta0_minus = beta0; beta0_minus(i) = beta0_minus(i) - delta;
            beta0_plus = beta0; beta0_plus(i) = beta0_plus(i) + delta;
            dfun_beta0(:,i) = real(fun(mu0,nu0,beta0_plus) - ...
                fun(mu0,nu0,beta0_minus))/2/delta;
        end
    end
else
    dfun_beta0  = funDiff_beta(mu0,nu0,beta0);
end

% B1 = \partial fun / \partial (mu,nu) = [funDmu0 funDnu0]
B10 = [dfun_mu0 dfun_nu0];

% B2 = \partial fun / \partial (beta) = funDbeta0
B20 = dfun_beta0;

% b  = fun(mu0,nu0,beta0)
b0  = fun_0;

% Jacobian matrix J = d(residual(x,y,mu,nu))/d(\mu,beta) =
% J = [-I 0; -d(fun(\mu,nu,beta))/d(mu0)  -d(fun(\mu,nu,beta))/d(beta0)]
J0  = [-diag(ones(n_mu0,1)) zeros(n_mu0,n_beta0); -dfun_mu0 -dfun_beta0];
end