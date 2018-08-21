function [excess, results] = island_compute_soe(prices,params,flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the small open economy version of Lyon Waugh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This first checks to see if it is in "callibration mode" or not
check_n_prices = length(prices);

if check_n_prices == (length(params.trans_mat)+1) + 3

    %disp('calibration mode')
    
    p_guess = prices;
    lambda = prices(end-3);
    params.labor_disutility = prices(end-2);
    params.m = prices(end-1);
    params.trade_cost = prices(end);
    prices = prices(1:end-4);

elseif check_n_prices == (length(params.trans_mat)+1) + 1
    p_guess = prices;
    lambda = prices(end-1);
    params.trade_cost = prices(end);
    prices = prices(1:end-2);
    
elseif check_n_prices == (length(params.trans_mat)+1)
    p_guess = prices;
    lambda = prices(end); 
    prices = prices(1:end-1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trans_mat = params.trans_mat;

grid = params.grid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalize aggregate price index to the value one and then compute pre and
% post tax wages...

ces_price = sum(params.invar_trans_mat.*prices.^(1-params.theta)).^(1./(1-params.theta));

prices = prices./ces_price;

params.world_prices = params.world_prices./ces_price;

pre_tax_wages = (prices.*params.shocks(:,1));

post_tax_wages = lambda.*(pre_tax_wages).^(1-params.tp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sovles the households problem...

[policy_assets, policy_move, value_fun] = mobility_labor_supply_smth(params,post_tax_wages,params.bsmooth);
% This is matching up with Spencers code...

% Compute the invariant distribution (across assets and locations) given
% the solution to the workers problem...

invariant_distribution = island_invariant(policy_assets,policy_move,trans_mat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Suming down a column gives the mass of workers in that location, or labor
% supply...
moving = sum(invariant_distribution.*policy_move(:,:,2) + invariant_distribution.*policy_move(:,:,4))';

labor_supply = sum(invariant_distribution.*policy_move(:,:,1) + invariant_distribution.*policy_move(:,:,2))';

locations = sum(invariant_distribution)';

compute_consumption
% Need to fix this, not a clear part of the code...

big_c = avg_consumption + agg_taxes;
% Update: With tariffs, I think this is the problem and equation (25). The
% deal is this is as if all government consumption comes from domestic
% production? But some of government consumption comes from tariff revenue
% which is external. 


[excess, stats] = market_clearing_soe(prices,labor_supply,big_c,agg_assets,moving,params, flag);

stats.moving = sum(moving);

gov_spending = agg_taxes + stats.tariff_revenue;

productivity = stats.output_hour;
output = stats.output;

avg_tax = sum(pre_tax_wages.*labor_supply)./sum(post_tax_wages.*labor_supply);

prob_assets = sum(invariant_distribution,2);

stats.frac_zero = sum(prob_assets(params.asset_space <= 0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
not_bounded = 10.*ones(check_n_prices - (length(params.trans_mat)),1);

total_import_trade_cost = (params.import_tau_shift.*params.trade_cost.*(1 + params.tariff));
%this is the total import trade costs, the shifted iceberg and then the
%tariff.

upper_bound = [total_import_trade_cost.*params.world_prices; not_bounded.^10 ]; 
% This is about the import side.

lower_bound = [params.world_prices./params.trade_cost; not_bounded.^-10];
% This is about the export side. 

%excess = min( max(-excess, log(p_guess) - log(upper_bound)), log(p_guess) - log(lower_bound));
excess(end+1) = agg_taxes./stats.output - params.gov;


if check_n_prices == (length(params.trans_mat)+1) + 3

    excess(end+1) = log(stats.labor_supply) - log(params.labor_supply);
    excess(end+1) = log(stats.moving) - log(params.mobility);
    excess(end+1) = log(stats.import_share) - log(params.trade_share);
    
elseif check_n_prices == (length(params.trans_mat)+1) + 1
    
    excess(end+1) = log(stats.import_share) - log(params.trade_share);
    
end
% Target this as a fraction of output

first_phi = fischer_burmeister((upper_bound) - (p_guess), excess);
excess = fischer_burmeister((p_guess) - (lower_bound), first_phi);
    
trade = stats.import_share;
    
ls = stats.labor_supply;
    
mov = sum(moving);

workers = labor_supply./params.invar_trans_mat;

avg_worker = sum(workers.*params.invar_trans_mat);
avg_wages = sum(pre_tax_wages.*params.invar_trans_mat);

OPterm1 =  avg_worker.*avg_wages;
OPterm2 = sum((avg_wages - pre_tax_wages).*(avg_worker - workers).*params.invar_trans_mat);


log_avg_wages = sum(log(pre_tax_wages).*params.invar_trans_mat);
var_wages = sum((log_avg_wages - log(pre_tax_wages)).^2.*params.invar_trans_mat);

marginal_rates = 1 - (1-params.tp).*lambda.*(pre_tax_wages).^(-params.tp);

[~, idx] = sort(marginal_rates);

income_prct = cumsum(params.invar_trans_mat(idx));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results.prices = prices;

results.trade = trade;

results.ls = ls;

results.mov = mov;

results.output = output;

results.opterm = OPterm2;

results.marginal_rates = marginal_rates(idx);

results.income_prct = income_prct;

results.welfare = welfare;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This plots and reports stuff. Is only called if flag equals two...
if flag ~=0

    [~, stats] = market_clearing_soe(prices,labor_supply,big_c,agg_assets,moving,params, flag);
    
    stats.moving = mov;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ces_price = sum(params.invar_trans_mat.*prices.^(1-params.theta)).^(1./(1-params.theta));
    
    prob_assets = sum(invariant_distribution,2);

    stats.frac_zero = sum(prob_assets(params.asset_space <= 0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('Mean Absolute Error in Excess Demand')
    disp(max(abs(excess)))
    disp('')
    disp('')
    disp('Real Interest Rate')
    disp([params.R, ces_price])
    disp('Average Tax Rate')
    disp(avg_tax)
    disp('Labor Supply, Mobility, Import Share, Export Share, Government Spending')
    disp([stats.labor_supply, stats.moving, stats.import_share, stats.export_share, (gov_spending./stats.output), agg_taxes./stats.output] )
    disp('Output, Output Per Hour')
    disp([stats.output, stats.output_hour] )
    disp('Welfare, Consumption')
    disp([welfare, avg_consumption, ] )
    disp('Trade Imbalance: Three ways')
    disp([stats.trade_balance_one, stats.trade_balance_two, stats.trade_balance_three] )
    disp('Aggregate Asset Holdings, Fraction Less Than Zero')
    disp([agg_assets, stats.frac_zero])
  
    disp('Output, First Term, Covraiance Term')
    disp([sum(pre_tax_wages.*labor_supply), OPterm1, OPterm2, var_wages ])
    %disp([corr(log(population(:)), log(params.shocks(:))), var(log(pre_tax_wages(:)))])
    
    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if flag == 1

    close all

    population  = reshape(locations,params.nshocks,params.nshocks);

    surf((population.*params.nshocks),'CDataMapping','scaled')
    set(gca,'XTick',1:1:params.nshocks)
    set(gca,'xticklabels',round(unique(params.shocks),2),'fontweight','bold')
    xlabel('Home Productivity')
    set(gca,'YTick',1:1:params.nshocks)
    set(gca,'yticklabels',round(unique(params.world_prices),2),'fontweight','bold')
    ylabel('World Prices')
    zlabel('Population (Normalized)')
    zlim([0,1])
    map = brewermap(10,'Blues');
    colormap(map)
    title('Distribution of Population Across Islands')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    ls_pop  = reshape(labor_supply./locations,params.nshocks,params.nshocks);
    surf((ls_pop),'CDataMapping','scaled')
    set(gca,'XTick',1:1:params.nshocks)
    set(gca,'xticklabels',round(unique(params.shocks),2),'fontweight','bold')
    xlabel('Home Productivity')
    set(gca,'YTick',1:1:params.nshocks)
    set(gca,'yticklabels',round(unique(params.world_prices),2),'fontweight','bold')
    ylabel('World Prices')
    zlabel('Fraction Working')
    zlim([0,1.5])
    map = brewermap(10,'Blues');
    colormap(map)
    title('Labor Supply')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    cons_dist  = reshape(consumption_by_location,params.nshocks,params.nshocks);

    surf(log(cons_dist))
    set(gca,'XTick',1:1:params.nshocks)
    set(gca,'xticklabels',round(unique(params.shocks),2),'fontweight','bold')
    xlabel('Home Productivity')
    set(gca,'YTick',1:1:params.nshocks)
    set(gca,'yticklabels',round(unique(params.world_prices),2),'fontweight','bold')
    %zlim([0,3])
    ylabel('World Prices')
    zlabel('Log Consumption')
    map = brewermap(10,'Blues');
    colormap(map)
    title('Consumption')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    ast_dist  = reshape(assets_by_location,params.nshocks,params.nshocks);

    surf((ast_dist))
    set(gca,'XTick',1:1:params.nshocks)
    set(gca,'xticklabels',round(unique(params.shocks),2),'fontweight','bold')
    xlabel('Home Productivity')
    set(gca,'YTick',1:1:params.nshocks)
    set(gca,'yticklabels',round(unique(params.world_prices),2),'fontweight','bold')
    zlim([params.grid(2),params.grid(3)])
    ylabel('World Prices')
    zlabel('Assets')
    map = brewermap(10,'Blues');
    colormap(map)
    title('Asset Holdings')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    wage_dist  = reshape(pre_tax_wages,params.nshocks,params.nshocks);


    surf(log(wage_dist))
    set(gca,'XTick',1:1:params.nshocks)
    set(gca,'xticklabels',round(unique(params.shocks),2),'fontweight','bold')
    xlabel('Home Productivity')
    set(gca,'YTick',1:1:params.nshocks)
    set(gca,'yticklabels',round(unique(params.world_prices),2),'fontweight','bold')
    %zlim([0,3])
    ylabel('World Prices')
    zlabel('Log (Pre-Tax) Wages')
    map = brewermap(10,'Blues');
    colormap(map)
    title('Wages')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    homeshare_dist  = reshape(stats.home_share,params.nshocks,params.nshocks);

    surf((homeshare_dist))
    set(gca,'XTick',1:1:params.nshocks)
    set(gca,'xticklabels',round(unique(params.shocks),2),'fontweight','bold')
    xlabel('Home Productivity')
    set(gca,'YTick',1:1:params.nshocks)
    set(gca,'yticklabels',round(unique(params.world_prices),2),'fontweight','bold')
    %zlim([0,2])
    ylabel('World Prices')
    zlabel('Home Share ')
    map = brewermap(10,'Blues');
    colormap(map)
    title('Home Share')


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure_density = figure;
    asset_space = params.asset_space;

    axes1 = axes('Parent',figure_density,'YGrid','on','XGrid','on','FontWeight','bold',...
        'FontSize',14);
    xlim([min(asset_space)-.05,max(asset_space)+.05]);

    hold(axes1,'all');

    xlabel('Assets','FontWeight','bold','FontSize',16);
    ylabel('Probability','FontWeight','bold','FontSize',16);
    prob_assets = sum(invariant_distribution,2);
    plot(asset_space,prob_assets,'LineWidth',3,'LineStyle','-','Color',[1 0 0])
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     home_shocks = params.shocks(:,1);
%     
%     stats.home_share = stats.home_share.^(params.theta);
%     
%     adh_data_all = [pre_tax_wages(params.sim_shocks), stats.home_share(params.sim_shocks),...
%         stats.per_worker(params.sim_shocks), home_shocks(params.sim_shocks), stats.imports(params.sim_shocks)];
%     
%     n_czs = 1000;
%     time_diff = 10;
% 
%     cm_zone = randperm(length(adh_data_all)-time_diff,n_czs)';
%     
%     adh_data_y1 = adh_data_all(cm_zone,:);
%     adh_data_y10 = adh_data_all(cm_zone + time_diff,:);
%     
%     year1 = 1990.*ones(length(cm_zone),1);
%     year2 = 2000.*ones(length(cm_zone),1);
%     
%     adh_data = [cm_zone, adh_data_y1, year1; cm_zone, adh_data_y10, year2];
%     
%     %csvwrite('C:\Users\mwaugh.NYC-STERN\Documents\GitHub\tradeexposure\model\matlab\simulated_data\adh_sim_data.csv',adh_data)
%     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     [b_level] = regress((log(pre_tax_wages(params.sim_shocks))), [ones(length(params.sim_shocks),1),...
%         (log(stats.home_share(params.sim_shocks))),(log(stats.per_worker(params.sim_shocks))), log(home_shocks(params.sim_shocks))]);
%     
%     [b_diff] = regress(diff(log(pre_tax_wages(params.sim_shocks))), [ones(length(params.sim_shocks)-1,1), diff(log(stats.home_share(params.sim_shocks))),...
%        diff(log(home_shocks(params.sim_shocks))), diff(log(stats.per_worker(params.sim_shocks)))]);
%    
%     [b_diff_adh] = regress(diff(log(pre_tax_wages(params.sim_shocks)),10), [ones(length(params.sim_shocks)-10,1),...
%         diff((stats.imports(params.sim_shocks)),10)./stats.per_worker(params.sim_shocks(1:end-10))]);
%     
% 
%     
%     b2sls_one = regress(diff(log(stats.home_share(params.sim_shocks))), [ones(length(params.sim_shocks)-1,1), diff(log(params.world_prices(params.sim_shocks)))]);
%     
%     x1hat = [ones(length(params.sim_shocks)-1,1), diff(log(params.world_prices(params.sim_shocks)))]*b2sls_one;
%           
%     [b_level] = regress(diff(log(pre_tax_wages(params.sim_shocks))), [ones(length(params.sim_shocks)-1,1), x1hat]);
%     % Issue seems to be that world prices and the home shocks are not
%     % orthogonal??? Not clear why that is?
%     
% 
%     
%     test_sim = 1;
    
    
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

