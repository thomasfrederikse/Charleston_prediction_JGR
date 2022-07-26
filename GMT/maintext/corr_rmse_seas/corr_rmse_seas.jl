# ACC and RMSE per season
using NCDatasets
using Statistics
using DelimitedFiles

fh = Dataset(homedir()*"/Data/ECCO/RISE/results/Charleston/Charleston_ccsm_TG_stats.nc","r")

rmse_pred = 100 .* fh["rmse_conv_obs_pred_mean_init_seas"][:]
rmse_seas = 100 .* fh["rmse_conv_obs_seas_mean_init_seas"][:]
rmse_pers = 100 .* fh["rmse_persistence_mean_init_seas"][:]
rmse_ssh = 100 .* fh["rmse_sla_pred_mean_init_seas"][:]
rmse_mean = 100 .* fh["rmse_obs_seas_mean_init_seas"][:]

corr_pred = fh["corr_conv_obs_pred_mean_init_seas"][:]
corr_seas = fh["corr_conv_obs_seas_mean_init_seas"][:]
corr_pers = fh["corr_persistence_mean_init_seas"][:]
corr_ssh = fh["corr_sla_pred_mean_init_seas"][:]

seasons = corr_persistence = fh["seasons"][:]

dir_gmt = homedir()*"/Projects/2020_ECCO_adjoint/GMT/Charleston/maintext/corr_rmse_seas/"
tlag = [1:12...]
for seas in 1:length(seasons)
    writedlm(dir_gmt*"rmse_pred_"*seasons[seas]*".txt",[tlag rmse_pred[seas,:]])
    writedlm(dir_gmt*"rmse_seas_"*seasons[seas]*".txt",[tlag rmse_seas[seas,:]])
    writedlm(dir_gmt*"rmse_pers_"*seasons[seas]*".txt",[tlag rmse_pers[seas,:]])
    writedlm(dir_gmt*"rmse_ssh_" *seasons[seas]*".txt",[tlag rmse_ssh[seas,:] ])
    writedlm(dir_gmt*"rmse_mean_"*seasons[seas]*".txt",[tlag rmse_mean[seas,:]])

    writedlm(dir_gmt*"corr_pred_"*seasons[seas]*".txt",[tlag corr_pred[seas,:]])
    writedlm(dir_gmt*"corr_seas_"*seasons[seas]*".txt",[tlag corr_seas[seas,:]])
    writedlm(dir_gmt*"corr_pers_"*seasons[seas]*".txt",[tlag corr_pers[seas,:]])
    writedlm(dir_gmt*"corr_ssh_"*seasons[seas]*".txt",[tlag corr_ssh[seas,:]])
end
close(fh)


