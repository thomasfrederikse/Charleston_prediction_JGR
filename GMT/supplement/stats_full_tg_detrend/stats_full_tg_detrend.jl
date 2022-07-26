# stats_full tide gauge
using DelimitedFiles
using NetCDF

fn_tg  = homedir() * "/Data/ECCO/RISE/results/Charleston/Charleston_ccsm_TG_detrend_stats.nc"
dir_gmt = homedir() * "/Projects/2020_ECCO_adjoint/GMT/Charleston/supplement/stats_full_tg_detrend/"

rmse_mss = (ncread(fn_tg,"rmse_conv_seas_mean")) .* 100
rmse_ecco_obs = (ncread(fn_tg,"rmse_conv_obs_mean")) .* 100
rmse_ecco_obs_seas = (ncread(fn_tg,"rmse_conv_obs_seas_mean")) .* 100
rmse_ccsm = (ncread(fn_tg,"rmse_conv_obs_pred_ens")) .* 100
rmse_ccsm_mn = (ncread(fn_tg,"rmse_conv_obs_pred_mean")) .* 100
rmse_persistence = (ncread(fn_tg,"rmse_persistence_mean")) .* 100
rmse_ccsm_sla = (ncread(fn_tg,"rmse_sla_pred_mean")) .* 100

corr_ecco_obs = (ncread(fn_tg,"corr_conv_obs_mean"))
corr_ecco_obs_seas = (ncread(fn_tg,"corr_conv_obs_seas_mean"))
corr_ccsm = (ncread(fn_tg,"corr_conv_obs_pred_ens"))
corr_ccsm_mn = (ncread(fn_tg,"corr_conv_obs_pred_mean")) 
corr_persistence = (ncread(fn_tg,"corr_persistence_mean")) 
corr_ccsm_sla = (ncread(fn_tg,"corr_sla_pred_mean")) 

writedlm(dir_gmt*"rmse_ccsm_sla.txt",[[1:12...] rmse_ccsm_sla],';')
writedlm(dir_gmt*"rmse_mss.txt",[[1:12...] rmse_mss],';')
writedlm(dir_gmt*"rmse_ecco_obs.txt",[[1:12...] rmse_ecco_obs],';')
writedlm(dir_gmt*"rmse_ecco_obs_seas.txt",[[1:12...] rmse_ecco_obs_seas],';')
writedlm(dir_gmt*"rmse_ccsm_mn.txt",[[1:12...] rmse_ccsm_mn],';')
for i ∈ 1:10
    writedlm(dir_gmt*"rmse_ccsm_"*string(i)*".txt",[[1:12...] rmse_ccsm[:,i]],';')
end
writedlm(dir_gmt*"rmse_persistence.txt",[[1:12...] rmse_persistence],';')

writedlm(dir_gmt*"corr_ccsm_sla.txt",[[1:12...] corr_ccsm_sla],';')

writedlm(dir_gmt*"corr_ecco_obs.txt",[[1:12...] corr_ecco_obs],';')
writedlm(dir_gmt*"corr_ecco_obs_seas.txt",[[1:12...] corr_ecco_obs_seas],';')
writedlm(dir_gmt*"corr_ccsm_mn.txt",[[1:12...] corr_ccsm_mn],';')
for i in 1:10
    writedlm(dir_gmt*"corr_ccsm_"*string(i)*".txt",[[1:12...] corr_ccsm[:,i]],';')
end
writedlm(dir_gmt*"corr_persistence.txt",[[1:12...] corr_persistence],';')
