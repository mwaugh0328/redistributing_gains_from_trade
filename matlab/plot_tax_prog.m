function [results] = plot_tax_prog(trade_share)

tax1 = -.25:.05:0;
tax2 = 0:.01:0.40;
tax3 = 0.35:0.05:0.60;
tax = [tax1, tax2(2:end), tax3(2:end)];
wel = zeros(length(tax),8);

tau = island_solve_counterfact_trade_NAG(0.18,trade_share);

parfor xxx = 1:length(tax)
    
    [out_results]  = island_solve_progresive_NAG(tax(xxx),tau,0,1);
    
    wel(xxx,:) = [out_results.trade_costs, out_results.trade, out_results.ls,...
        out_results.mov, out_results.output, out_results.opterm, out_results.welfare, double(out_results.exit_flag)];

    disp(wel(xxx,:))
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = [tax', wel];

did_not_solve = (results(:,end) ~=0);

all_results = results(did_not_solve==0,:);

Y = all_results(:,8); X = [all_results(:,1)];

[bhat,~,~,~,stats] = regress(Y,[ones(length(X),1), X, X.^2]);

disp(stats)

X = [results(:,1)];

ypred = bhat'*[ones(length(X),1), X, X.^2]';

results = [tax', wel, ypred'];

file_name = strcat('results',num2str(trade_share),'.mat');

cd('.\plot_model_data')

save(file_name, 'results')

cd('..\')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%