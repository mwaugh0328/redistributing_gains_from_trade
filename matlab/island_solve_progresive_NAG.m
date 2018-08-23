function [results] = island_solve_progresive_NAG(tp,trade_cost,tariff,flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

%ifail = x06aa(int64(24));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Calibrated values....

load calibration
% Given the ordering of the calibration (labor disutility, moving cost)
% this passes through the correct values.
params.labor_disutility = calibration_results(1); 

params.m = calibration_results(2); % Moving cost about 3 percent

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set other parameter values...
params.R = 1.02; % Gross real interest rate

params.beta = 0.95; % Discount factor

params.theta = 4.0; % Elasticity of substitution on the CES demand cure

params.home_production = 0.00; % Not sure what to do with this. Chang, etc. do not have this, set 0. What about government policy?

params.trade_cost = trade_cost; % About 13 percent of GDP.

params.tariff = tariff;

params.import_tau_shift = 1.00;

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

guess = exp(.25.*log(params.world_prices));
guess(end+1) = 0.80;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tic;  xxx = island_compute_soe((guess),params,0); toc
% disp(norm(xxx))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nval = length(guess);
diag_adjust = round(nval-1);

tic
[prices, fvec, diag, nfev,~,~,~,~,ifail] = c05qc(@fcn, log(guess), int64(diag_adjust),...
   int64(diag_adjust), int64(1), ones(nval,1), int64(5),'user', params,'epsfcn', 10^-5,'xtol',10^-5);
toc

% if flag == 0
% 
%    [~, out_results]= island_compute_soe(exp(prices),params,0);
% 
%    results = -out_results.welfare;
% 
% else

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('')
disp('')
disp('Number of function evaluations')
disp(nfev)
disp('Error Warning')
disp(ifail)

exit_flag = ifail;

%[excess, prices, trade, ls, mov, welfare]
[~, results]= island_compute_soe(exp(prices),params,flag);

results.exit_flag = exit_flag;
results.trade_costs = params.trade_cost;

%save solution prices


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fvec, user, iflag] = fcn(n, x, fvec, user, iflag)
 if iflag ~=0
   fvec = zeros(n, 1);
   fvec(1:n) = island_compute_soe(exp(x),user,0);
 else
   fprintf('objective = %10e\n', max(abs(fvec)))
 end
