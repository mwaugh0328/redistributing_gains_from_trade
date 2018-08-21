function invariant_distribution = island_invariant(policy_assets, policy_move,trans_mat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we want compute the invariance distribution over states and assets.
% This is the trickiest part I found, the approach here is to contruct the
% markov transition matrix for the (asset,state) combination. Then figure
% out the invariant distribution...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n_asset_states,n_shocks,~] = size(policy_assets);

QA_mat = sparse(repmat(1:n_asset_states*n_shocks,1,4),policy_assets(:)',policy_move(:)',n_asset_states*n_shocks,n_asset_states);

% Now what this is going to do is to construct a matrix which is
% (n_assets*n_shocks)*n_assets. Each entry will represent: given an
% asset,shock state the probability that I arrive at a particular asset state
% The sparse command is just doing this is a quick way... need to look in
% documentation to see exactly what this is doing.
invar_trans_mat = trans_mat^50000;
invar_trans_mat = invar_trans_mat(1,:);

no_move_trans = repmat(trans_mat,n_asset_states,1);

policy_move_probs = reshape(policy_move,n_shocks*n_asset_states,4);

%QW_mat = zeros(n_asset_states*n_shocks,n_shocks);

asset_state_transition = zeros(n_asset_states*n_shocks,n_asset_states*n_shocks);

shock_index = repmat(1:n_shocks,n_asset_states,1);
shock_index = shock_index(:);

for zzz = 1:n_shocks*n_asset_states

    prob_move = sum(policy_move_probs(zzz,[2,4]));
    prob_stay = sum(policy_move_probs(zzz,[1,3]));

    QW_mat = prob_stay.*trans_mat(shock_index(zzz),:) + prob_move.*invar_trans_mat;

    asset_state_transition(zzz,:) = kron(QW_mat,QA_mat(zzz,:));

end

%[L] = invar_mplt(asset_state_transition);

%L = invar_mplt_mex(asset_state_transition);

% %L = 1./(n_shocks*n_asset_states).*ones(1,n_shocks*n_asset_states);
%
% % the strategy is to compute Lnew = Q' * L; until Lnew -  L.
% % This is similar to doing mpower, but will stop whenever we are "close
% % enough" and doesn't do quite as much computation.

asset_state_transition = sparse(asset_state_transition);

L = zeros(1,n_shocks*n_asset_states);
L((n_shocks*n_asset_states)/2) = 1;

for zzz = 1:2000
    L_new = L*asset_state_transition;
    %L_new = mtimesx(L,asset_state_transition,'SPEED');

    %[L_new, ~, ~, ~] = f01ck(L, asset_state_transition, opt);

    if norm(L_new-L) < 10^-10
        break
    end

    L = L_new;
end

%step_one_invariant = L;

invariant_distribution = reshape(L,n_asset_states,n_shocks);


% this is correct given the new ordering...
