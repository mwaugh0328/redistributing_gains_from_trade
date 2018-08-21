island_solve_calibrate_NAG(0.18)

trade_share = [0.05,0.10,0.20,0.30,0.40];

for xxx = 1:length(trade_share)
    
    plot_tax_prog(trade_share(xxx))
    
end
    