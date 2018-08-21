%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to generate the markov chain over the local productivity shocks and
% world prices....

n_shocks = params.nshocks;

[shocks_local,trans_mat_local] = rouwenhorst(n_shocks, p, std_inv);

[shocks_world,trans_mat_world] = rouwenhorst(n_shocks, p, std_inv_world);

params.trans_mat = kron(trans_mat_local, trans_mat_world);

[shock_home, world_price]=meshgrid(shocks_local,shocks_world');

forg_index = flipud(rot90(reshape((1:n_shocks^2),n_shocks,n_shocks)));
params.forg_index = forg_index(:);

params.shocks = exp(shock_home(:));
params.world_prices = exp(world_price(:));

invar_trans_mat = params.trans_mat^50000;
params.invar_trans_mat = invar_trans_mat(1,:)';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This tests the markov chain simmulation...

%rng(03281978)

%[~, params.sim_shocks] = hmmgenerate(10000,params.trans_mat,ones(length(params.trans_mat),1));

% rho = regress(log(params.shocks(shock_states_p(2:end))),log(params.shocks(shock_states_p(1:end-1))));
% 
% disp(rho)
% 
% test = 1;
