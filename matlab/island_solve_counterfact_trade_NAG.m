function [results] = island_solve_counterfact_trade_NAG(tp,tradeshare)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

%ifail = x06aa(int64(24));

load calibration

params.gamma = 2; % Risk aversion
params.R = 1.02; % Gross real interest rate
params.beta = 0.95; % Discount factor

params.labor_disutility = calibration_results(1); % Chang, etc. 60 percent employment ratio. Page 1943

params.m = calibration_results(2); % Moving cost about 3 percent

params.theta = 4; % Elasticity of substitution on the CES demand cure

params.home_production = 0.00; % Not sure what to do with this. Chang, etc. do not have this, set 0. What about government policy?

params.trade_share = tradeshare; % About 10 percent of GDP in 1990 % Same deal.

params.tariff = 0.00;

params.import_tau_shift = 1.00;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The tax function... disposable income = lambda^(1-tp)

params.gov = 0.19;

params.tp = tp;

params.bsmooth = 0.20; % The smoothing parameter.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shocks and the grid...
params.nshocks = 10;
params.grid = [50, -0.45, 6];
params.asset_space =  clustergrid(params.grid(2), 0.05, params.grid(3), 20, 30, 0.7, 0.3);
%params.asset_space =  linspace(params.grid(2),params.grid(3),params.grid(1));%

p = 0.95;
std_inv = sqrt(0.017).*(params.theta/(params.theta - 1)); % Need to scale up as wages are dapened by theta...
% Now Kaplan (2012) and Ventura et.
% al report yearly values of 0.95 and sqrt(0.017). Second issue is that
% this is the wage process, the world price process.

rel_var = 1.0; 
% This adjusts the relative variance of the world prices.
std_inv_world = rel_var.*std_inv;

open_economy_markovchain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load solution
%guess = exp(prices);
%
guess = exp(0.30.*log(params.world_prices));
guess(end+1) = 1.0;
guess(end+1) = 3.0; % This is the trade costs

%guess = 0.95.*guess + 0.05.*guess_hat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tic;  xxx = island_compute_soe((guess),params,0); toc
% disp(norm(xxx))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nval = length(guess);
diag_adjust = round(nval-3);

tic
[prices, fvec, diag, nfev,~,~,~,~,ifail] = c05qc(@fcn, log(guess), int64(diag_adjust),...
   int64(diag_adjust), int64(1), ones(nval,1), int64(5),'user', params,'epsfcn', 10^-4,'xtol',10^-5);
toc

% options = optimoptions('fsolve', 'Display','iter','MaxFunEvals',3000,'MaxIter',...
%     1e2,'TolX',1e-5,'UseParallel',true,'Algorithm','trust-region','FiniteDifferenceType','central')
% % 'SubproblemAlgorithm','cg'
% % options = optimoptions('fsolve', 'Display','iter','MaxFunEvals',3000,'MaxIter',...
% %     1e2,'TolFun',1e-2,'TolX',1e-2,'UseParallel',true,'Algorithm','levenberg-marquardt','InitDamping',0.001)
% 
% tic
% [prices, ~, exit_flag] = fsolve(@(xxx) island_compute_soe(exp(xxx),params,0),log(guess),options);
% toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('')
disp('')
disp('Number of function evaluations')
disp(nfev)
disp('Error Warning')
disp(ifail)

exit_flag = ifail;

% [~, results]= island_compute_soe(exp(prices),params,0);
% disp(trade)

results = exp(prices(end));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fvec, user, iflag] = fcn(n, x, fvec, user, iflag)
 if iflag ~=0
   fvec = zeros(n, 1);
   fvec(1:n) = island_compute_soe(exp(x),user,0);
 else
   fprintf('objective = %10e\n', max(abs(fvec)))
 end
