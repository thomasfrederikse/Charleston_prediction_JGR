# Compare altimetry and tide gauges
using DelimitedFiles
using NetCDF
using Plots
using Statistics


function main()
    fn_alt  = homedir() * "/Data/ECCO/RISE/results/Charleston/Charleston_Alt_tseries.nc"
    fn_tg  = homedir() * "/Data/ECCO/RISE/results/Charleston/Charleston_Alt_tseries.nc"
    dir_gmt = homedir() * "/Projects/2020_ECCO_adjoint/GMT/supplement/alt_tg_comp/"

    time = ncread(fn_alt,"time")[1,:] ./ 365.24 .+ 1900

    alt = (ncread(fn_alt,"obs_Alt")[1,:]) .* 100
    tg = (ncread(fn_tg,"obs_TG")[1,:]) .* 100

    alt_seas = (ncread(fn_alt,"obs_Alt_seas")[1,:]) .* 100
    tg_seas = (ncread(fn_tg,"obs_TG_seas")[1,:]) .* 100

    alt_noseas = alt - alt_seas
    tg_noseas = tg - tg_seas

    alt = alt .- mean(alt)
    tg = tg .- mean(tg)
    alt_noseas = alt_noseas .- mean(alt_noseas)
    tg_noseas = tg_noseas .- mean(tg_noseas)
    
    writedlm(dir_gmt*"alt.txt",[t alt],';')
    writedlm(dir_gmt*"tg.txt",[t tg],';')
    writedlm(dir_gmt*"alt_noseas.txt",[t alt_noseas],';')
    writedlm(dir_gmt*"tg_noseas.txt",[t tg_noseas],';')

    writedlm(dir_gmt*"tg_vs_alt_noseas.txt",[detrend(alt_noseas) detrend(tg_noseas)],';')
    writedlm(dir_gmt*"tg_vs_alt.txt",[detrend(alt) detrend(tg)],';')

    writedlm(dir_gmt*"straightline.txt",[[-30:30...] [-30:30...]],';')

    plot(t,alt,label="Altimetry",linewidth=2,legend=:topleft)
    plot!(t,tg,label="Tide gauge",linewidth=2)
    ylabel!("Height (mm)")
    savefig("alt_tg_ts.png")

    plot(t,alt_noseas,label="Altimetry",linewidth=2,legend=:topleft)
    plot!(t,tg_noseas,label="Tide gauge",linewidth=2)
    ylabel!("Height (mm)")
    savefig("alt_tg_ts_noseas.png")

    plot(detrend(alt),detrend(tg),marker=:circ,linewidth=0,label="")
    plot!([-20:30...],[-20:30...],label="y = x",color=:black,legend=:topleft)
    title!("Seasonal cycle retained")
    xlabel!("Altimetry (mm)")
    ylabel!("Tide gauge (mm)")
    savefig("alt_tg_hist.png")

    plot(detrend(alt_noseas),detrend(tg_noseas),marker=:circ,linewidth=0,label="")
    title!("Seasonal cycle removed")
    plot!([-15:20...],[-15:20...],label="y = x",color=:black,legend=:topleft)
    xlabel!("Altimetry (mm)")
    ylabel!("Tide gauge (mm)")
    savefig("alt_tg_hist_noseas.png")

    cor(detrend(alt),detrend(tg))
    cor(detrend(alt_noseas),detrend(tg_noseas))
end

function detrend(ts)
    # Design matrix
    amat = ones(size(ts,1),2);
    amat[:,2] = [1:length(ts)...];
    amat_tr = transpose(amat)
    amat_sq = inv(amat_tr*amat)
    sol = amat_sq * (amat_tr*ts)
    ts_detrend = ts - amat*sol
    return ts_detrend
end

function lsqtrend(x,y)
    amat = ones(size(x,1),2);        
    amat[:,2] = x
    amat_tr = transpose(amat)
    amat_sq = inv(amat_tr*amat)
    trend = (amat_sq * (amat_tr*y))[2]
    return trend
end

