function [yyy] = plot_opt_mag(xxx)

tau_tax = [5.618187, 0.18; ...
2.332627, 0.18; ...    
2.332627, 0.27; ...
1.964557, 0.32; ...
1.749793, 0.37; ...
1.586763, 0.45;];

marg_rates = zeros(100,length(tau_tax));
incom_prct = zeros(100,length(tau_tax));

parfor xxx = 1:length(tau_tax)
    
    [out_results]  = island_solve_progresive_NAG(tau_tax(xxx,2),tau_tax(xxx,1),0,1);
    
    marg_rates(:,xxx) = [out_results.marginal_rates];
    incom_prct(:,xxx) = [out_results.income_prct];
    
end

cd('.\plot_model_data')

save('opt_marg_rates.mat', 'marg_rates')
save('opt_incom_prct.mat', 'incom_prct')

cd('..\')