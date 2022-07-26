# rmse_proc
using DelimitedFiles
using NetCDF

fn_tg  = homedir() * "/Projects/2020_ECCO_adjoint/Results/Charleston_TG_stats.nc"
dir_gmt = homedir() * "/Projects/2020_ECCO_adjoint/GMT/rmse_proc/"

rmse_proc_mss = (ncread(fn_tg,"rmse_proc_mss")) .* 100
rmse_proc_ccsm = (ncread(fn_tg,"rmse_proc_ccsm")) .* 100
rmse_proc_ecco_obs_seas = (ncread(fn_tg,"rmse_proc_ecco_obs_seas")) .* 100
rmse_proc_ccsm_mn = (ncread(fn_tg,"rmse_proc_ccsm_mn")) .* 100
rmse_proc_ccsm_mn = (ncread(fn_tg,"rmse_proc_ccsm_mn")) .* 100

procs = ["empmr","qnet","tauu","tauv"]
for (idx,proc) ∈ enumerate(procs)
    writedlm(dir_gmt*"rmse_"*proc*"_seas.txt",[[1:12...] rmse_proc_mss[:,idx]],';')
    writedlm(dir_gmt*"rmse_"*proc*"_mean.txt",[[1:12...] rmse_proc_ccsm_mn[:,idx]],';')
    writedlm(dir_gmt*"rmse_"*proc*"_ecco_obs_seas.txt",[[1:12...] rmse_proc_ecco_obs_seas[:,idx]],';')
    for ens ∈ 1:10
        writedlm(dir_gmt*"rmse_"*proc*"_"*string(ens)*".txt",[[1:12...] rmse_proc_ccsm[:,idx,ens]],';')
    end
end
