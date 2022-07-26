# comp_obs_ecco
using DelimitedFiles
using NCDatasets
using Statistics
using Dates

fh_rec_tg  = Dataset(homedir() * "/Data/ECCO/RISE/results/Charleston/Charleston_ccsm_TG_tseries.nc","r")
fh_rec_alt  = Dataset(homedir() * "/Data/ECCO/RISE/results/Charleston/Charleston_ccsm_Alt_tseries.nc","r")
dir_gmt = homedir() * "/Projects/2020_ECCO_adjoint/GMT/comp_obs_ecco/"

# Time 
rec_time = fh_rec_tg["time"][1,:]
rec_time = convert.(Float32,[Year(rec_time[t]).value + (Month(rec_time[t]).value)/12f0 - 1f0/24f0 for t in 1:length(rec_time)])

# Observations
obs_alt  = 100 .* fh_rec_alt["obs_Alt"][1,:]
obs_alt_seas  = 100 .* fh_rec_alt["obs_Alt_seas"][1,:]

obs_tg   = 100 .* fh_rec_tg["obs_TG"][1,:]
obs_tg_seas   = 100 .* fh_rec_tg["obs_TG_seas"][1,:]

# ECCO
sshdyn_alt = 100 .* fh_rec_alt["pred_sla_obs_mean"][1,:]
ssh_rec_alt = 100 .* fh_rec_alt["pred_conv_obs_mean"][1,:]

sshdyn_tg = 100 .* fh_rec_tg["pred_sla_obs_mean"][1,:]
ssh_rec_tg = 100 .* fh_rec_tg["pred_conv_obs_mean"][1,:]

### Remove MSC from everything
# obs_alt .-= obs_alt_seas
obs_tg .-= obs_tg_seas
sshdyn_tg .-= obs_tg_seas
ssh_rec_tg .-= obs_tg_seas

obs_alt .-= obs_alt_seas
sshdyn_alt .-= obs_alt_seas
ssh_rec_alt .-= obs_alt_seas


# t = ncread(fn_obs,"time")[:] ./ 365.24 .+ 1900
# alt = (ncread(fn_obs,"obs_msl")[1,:] .- ncread(fn_alt,"obs_msl_seas")[1,:]) .* 100
# tg = (ncread(fn_obs,"obs_msl")[1,:] .- ncread(fn_tg,"obs_msl_seas")[1,:]) .* 100

# sshdyn = (ncread(fn_tg,"SSHDYN_msl")[1,:] .- ncread(fn_tg,"SSHDYN_msl_seas")[1,:]) .* 100
# ssh_rec = (ncread(fn_tg,"rec_ecco_obs")[1,:] .- ncread(fn_tg,"obs_msl_seas")[1,:]) .* 100

# Correlations
cor(obs_alt,sshdyn_alt)
cor(obs_alt,ssh_rec_alt)
cor(obs_tg,sshdyn_tg)
cor(obs_tg,ssh_rec_tg)

# rsq
compute_rsq(obs_tg,sshdyn_tg)
compute_rsq(obs_tg,ssh_rec_tg)

# Without 2004
acc_idx = floor.(rec_time) .!= 2004.0f0
cor(obs_tg[acc_idx],sshdyn_tg[acc_idx])
cor(obs_tg[acc_idx],ssh_rec_tg[acc_idx])


# Save
writedlm(dir_gmt*"alt.txt",[t alt],';')
writedlm(dir_gmt*"tg.txt",[t tg],';')
writedlm(dir_gmt*"ECCO_adj.txt",[t ssh_rec],';')
writedlm(dir_gmt*"ECCO_ssh.txt",[t sshdyn],';')

function compute_rsq(truth,model)
    tminmod = truth - model
    rsq = 1 - (sum((tminmod .- mean(tminmod)).^2))/(sum((truth .- mean(truth)).^2))
    return rsq
end
