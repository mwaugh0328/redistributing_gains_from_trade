function [excess_demand, stats] = market_clearing_soe(prices,labor_supply,consumption,assets,moving,params,flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes labor demand and excess demand function as specified by a
% standard CES demand function...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First, compute all the paymenst to labor. This will compute aggregate
% income.

home_shocks = params.shocks(:,1);

n_markets = length(home_shocks);
% Note the prices here should always be the price including the trade cost
% if imported....so consumer prices. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wages = prices.*home_shocks;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First compute production

home_supply = home_shocks.*labor_supply;

% mass of islands times productivity times number of workers on island. I'm
% not sure if the probability weight of the shock should show up here.
% Implicitly, is this not in the labor supply already.

%supply = zeros(n_markets,1);

supply = home_supply;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ces_price = sum(params.invar_trans_mat.*prices.^(1-params.theta)).^(1./(1-params.theta));

real_expenditure = consumption;
% Note that when we are doing this, we are implicity cleaing the final
% goods market. Its ok. just want to be explicit.

demand = (prices./ces_price).^(-params.theta).*real_expenditure;

demand = params.invar_trans_mat.*demand;

excess_demand = (log(demand) - log(supply));

imports =  max(demand - home_supply,0)./(1-params.tariff) ;
% See how this is set up in paper, the difference between demand and supply
% is the imports net of tariffs. So to get total imports, divide by the net
% tariff rate.

stats.tariff_revenue = sum(prices.*params.tariff.*imports);
% A question is about the pi thing... its aleady their know?


exports = -min(demand - home_supply,0);

imported = imports ~= 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stats.consumption = consumption;

stats.consumption2 = sum(prices.*demand);

stats.output = sum(wages.*labor_supply);

stats.output2 = sum(prices.*home_supply);

stats.trade_balance_one = stats.output - (stats.consumption + stats.tariff_revenue);

stats.trade_balance_two = ((sum(prices.*exports))...
    - sum(prices.*imports));

stats.trade_balance_three = -assets.*(params.R-1) + ...
    sum(moving).*params.m - stats.tariff_revenue;

stats.import_share = sum(prices.*imports)./stats.output;

stats.export_share = sum(prices.*exports)./stats.output;

stats.imports = prices.*imports; 

stats.labor_supply = sum(labor_supply);

stats.output_hour = stats.output./stats.labor_supply;

if flag ~= 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stats.labor_supply = sum(labor_supply);

stats.output_hour = stats.output./stats.labor_supply;

stats.per_worker = labor_supply./params.invar_trans_mat;

stats.home_share = (prices.*home_supply./(prices.*demand)).^(1./params.theta);

stats.wage_hat = stats.home_share.*stats.per_worker.^(-1./params.theta).*...
    home_shocks.^((params.theta-1)./params.theta).*consumption.^(1/params.theta);
% This should be the same as wages (still a bit puzzled that I feel like I
% have an extra pi(z) floating around....
%     ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%











