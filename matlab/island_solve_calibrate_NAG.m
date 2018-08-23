function [results] = island_solve_calibrate_NAG(tp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

%ifail = x06aa(int64(24));

params.R = 1.02; % Gross real interest rate

params.beta = 0.95; % Discount factor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Values that we calibrate too...

params.labor_supply = 0.67; % This is the calibration routine, so we will target 67 percent participation while solving equillibrium

params.mobility = 0.03; % This is the calibration routine, so we will targiet mobility of 3 percent

params.trade_share = 0.10; % This is the calibration routine, About 10 percent of GDP in 1990 % Same deal.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.theta = 4.0; % Elasticity of substitution on the CES demand cure

params.home_production = 0.00; % Home production, for this paper we set to zero.

params.tariff = 0.00; % tariff if we put it in there. 

params.import_tau_shift = 1.00; % reltive trade cocst

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The tax function... disposable income = lambda^(1-tp)

params.gov = 0.19; % Government spending as a fraction of income.

params.tp = tp;

params.bsmooth = 0.20; % The smoothing parameter.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shocks and the grid...
params.nshocks = 10; % Set to 10 to exactly replicate paper.
params.grid = [50, -0.45, 6]; % This is hand picked to get 40 percent with < 0 assets.
params.asset_space =  clustergrid(params.grid(2), 0.05, params.grid(3), 20, 30, 0.7, 0.3);
%params.asset_space =  linspace(params.grid(2),params.grid(3),params.grid(1));%

p = 0.95;
std_inv = sqrt(0.017).*(params.theta/(params.theta - 1)); % Need to scale up as wages are dapened by theta...
% Now Kaplan (2012) and Ventura et.
% al report yearly values of 0.95 and sqrt(0.017). .

rel_var = 1.0; 
% This adjusts the relative variance of the world prices.
std_inv_world = rel_var.*std_inv;

% References another .m file to construct the markov change. 
open_economy_markovchain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

guess = exp(0.5.*log(params.world_prices));
guess(end+1) = 1.00;
guess(end+1) = 1.24; % This is the labor supply disutility
guess(end+1) = 1.0; % this is the migration cost...
guess(end+1) = 2.5; % This is the trade costs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Just to see what is initially going on...

tic;  xxx = island_compute_soe((guess),params,0); toc
disp(norm(xxx))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is perferred approach, use NAG solver...

nval = length(guess);
diag_adjust = round(nval-1);

tic
[prices, fvec, diag, nfev,~,~,~,~,ifail] = c05qc(@fcn, log(guess), int64(diag_adjust),...
   int64(diag_adjust), int64(1), ones(nval,1), int64(5),'user', params,'epsfcn', 10^-5,'xtol',10^-5);
toc

disp('')
disp('')
disp('Number of function evaluations')
disp(nfev)
disp('Error Warning')
disp(ifail)

%exit_flag = ifail;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is using Matlab's solver...

% options = optimoptions('fsolve', 'Display','iter','MaxFunEvals',6000,'MaxIter',...
%     1e2,'TolX',1e-6,'UseParallel',true,'Algorithm','trust-region-reflective','FiniteDifferenceType','central');
% 
% tic
% [prices, ~, exit_flag] = fsolve(@(xxx) island_compute_soe(exp(xxx),params,0),log(guess),options);
% toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate figures

%[excess, prices, trade, ls, mov, welfare]
[~, results]= island_compute_soe(exp(prices),params,1);

calibration_results = exp([prices(end-2:end)])';
disp('Calibration Results: Labor Supply Disutility; Moving Cost, Trade Costs')
disp(calibration_results)

save calibration calibration_results 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fvec, user, iflag] = fcn(n, x, fvec, user, iflag)
 if iflag ~=0
   fvec = zeros(n, 1);
   fvec(1:n) = island_compute_soe(exp(x),user,0);
 else
   fprintf('objective = %10e\n', max(abs(fvec)))
 end
