%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This computes consumption by location....

m = params.m;

grid = params.grid;

home_production = params.home_production;

R = params.R;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up grid...

n_asset_states = grid(1);
asset_space = params.asset_space;

consumption_option = zeros(size(policy_assets));
assets_option = zeros(size(policy_assets));
taxes_option = zeros(size(policy_assets));
labor_supply_option = zeros(size(policy_assets));


[~,~,n_options] = size(policy_move);

consumption = zeros(n_asset_states,length(post_tax_wages));
assets = zeros(n_asset_states,length(post_tax_wages));
taxes = zeros(n_asset_states,length(post_tax_wages));
welfare_function = zeros(n_asset_states,length(post_tax_wages));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% So first walk through this, wage state, asset state...
for yyy = 1:n_options
            
            move = 0;
            not_working = 0;
            
            if yyy == 3 || yyy == 4
                not_working = 1;
            end
            
            if yyy == 2 || yyy == 4
                move = 1;
            end

    for zzz = 1:length(post_tax_wages)
    
        for xxx = 1:n_asset_states
                          
        % Then given their asset state they get returns on that asset minus
        % their asset policy + labor earnings if they work, home production
        % if they don't net of moving if they do move....
        assets_option(xxx,zzz,yyy) = asset_space(policy_assets(xxx,zzz,yyy));
        
        consumption_option(xxx,zzz,yyy) = (R.*asset_space(xxx)- assets_option(xxx,zzz,yyy)...
                            + post_tax_wages(zzz).*(~not_working) ...
                            - m.*move);
                        
        taxes_option(xxx,zzz,yyy) = (pre_tax_wages(zzz) - post_tax_wages(zzz)).*(~not_working);
        
        labor_supply_option(xxx,zzz,yyy) = (~not_working);
                        
        end
    end
    
    consumption = consumption_option(:,:,yyy).*policy_move(:,:,yyy) + consumption;
    
    taxes = taxes_option(:,:,yyy).*policy_move(:,:,yyy) + taxes;
    
    assets = assets_option(:,:,yyy).*policy_move(:,:,yyy) + assets;
    
    welfare_function = value_fun(:,:,yyy).*policy_move(:,:,yyy)  + welfare_function;
    
end

% Then multiply the consumption matrix by the invariant distribution, then
% sum down the column (normalized by the number of guys in the column).
% This gives the ``average consumption'' in that state...

consumption_by_location = (sum(consumption.*invariant_distribution)./sum(invariant_distribution))';

avg_consumption = sum(sum(consumption.*invariant_distribution));

agg_taxes = sum(sum(taxes.*invariant_distribution));

assets_by_location = (sum(assets.*invariant_distribution)./sum(invariant_distribution))'; 

agg_assets = sum(sum(assets.*invariant_distribution));

welfare = sum(sum(welfare_function.*invariant_distribution));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%