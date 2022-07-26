# Example time series
using DelimitedFiles
using NCDatasets
using Dates

fn_tg  = homedir() * "/Data/ECCO/RISE/results/Charleston/Charleston_ccsm_Alt_tseries.nc"
dir_gmt = homedir() * "/Projects/2020_ECCO_adjoint/GMT/Charleston/supplement/example_rec_tseries_alt/"
mnum=1

fh = Dataset(fn_tg,"r")
tval = variable(fh,"time")[:] ./ 365.25 .+ 1900
tstart = variable(fh,"tstart")[:] ./ 365.25 .+ 1900
obs_seas = fh["pred_conv_obs_seas_mean"][:] .* 100
obs_pred = fh["pred_conv_obs_pred_mean"][:] .* 100
obs_pred_ens = fh["pred_conv_obs_pred_ens"][:] .* 100
sla_pred = fh["pred_sla_pred_mean"][:] .* 100
sla_pred_ens = fh["pred_sla_pred_ens"][:] .* 100
obs = fh["obs_TG"][:] .* 100
msc = fh["obs_TG_seas"][:] .* 100

tval_1m = tval[:,1:12:end][:]
obs_seas_1m = obs_seas[:,1:12:end][:]
obs_pred_1m = obs_pred[:,1:12:end][:]
obs_pred_ens_1m = zeros(Float32,10,264)
sla_pred_1m = sla_pred[:,1:12:end][:]
sla_pred_ens_1m = zeros(Float32,10,264)
for i in 1:10
    obs_pred_ens_1m[i,:] = obs_pred_ens[:,i,1:12:end][:]
    sla_pred_ens_1m[i,:] = sla_pred_ens[:,i,1:12:end][:]
end
obs_1m = obs[:,1:12:end][:]
msc_1m = msc[:,1:12:end][:]

tval_4m = tval[:,4:12:end][:]
obs_seas_4m = obs_seas[:,4:12:end][:]
obs_pred_4m = obs_pred[:,4:12:end][:]
obs_pred_ens_4m = zeros(Float32,10,264)
sla_pred_4m = sla_pred[:,4:12:end][:]
sla_pred_ens_4m = zeros(Float32,10,264)

for i in 1:10
    obs_pred_ens_4m[i,:] = obs_pred_ens[:,i,4:12:end][:]
    sla_pred_ens_4m[i,:] = sla_pred_ens[:,i,4:12:end][:]
end
obs_4m = obs[:,4:12:end][:]
msc_4m = msc[:,4:12:end][:]

# Save data
writedlm(dir_gmt*"obs_pred_1m.txt",[tval_1m obs_pred_1m],';')
writedlm(dir_gmt*"sla_pred_1m.txt",[tval_1m sla_pred_1m],';')
writedlm(dir_gmt*"obs_seas_1m.txt",[tval_1m obs_seas_1m],';')
writedlm(dir_gmt*"obs_1m.txt",[tval_1m obs_1m.-msc_1m],';')

writedlm(dir_gmt*"obs_pred_4m.txt",[tval_4m obs_pred_4m],';')
writedlm(dir_gmt*"sla_pred_4m.txt",[tval_4m sla_pred_4m],';')
writedlm(dir_gmt*"obs_seas_4m.txt",[tval_4m obs_seas_4m],';')
writedlm(dir_gmt*"obs_4m.txt",[tval_4m obs_4m.-msc_4m],';')
for i in 1:10
    writedlm(dir_gmt*"obs_pred_ens_"*string(i)*"_1m.txt",[tval_1m obs_pred_ens_1m[i,:]])
    writedlm(dir_gmt*"obs_pred_ens_"*string(i)*"_4m.txt",[tval_4m obs_pred_ens_4m[i,:]])
    writedlm(dir_gmt*"sla_pred_ens_"*string(i)*"_1m.txt",[tval_1m sla_pred_ens_1m[i,:]])
    writedlm(dir_gmt*"sla_pred_ens_"*string(i)*"_4m.txt",[tval_4m sla_pred_ens_4m[i,:]])
end

