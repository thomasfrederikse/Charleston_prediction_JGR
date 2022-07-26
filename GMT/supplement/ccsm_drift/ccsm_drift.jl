# ----------------------------------------------------------------
# Figure with 4 panels:
# 1. time series of CCSM4 and CCSM4 dedrifted ssh vs observations
# 2. Spatial pattern of drift in first year_sens
# 3 and 4. RMSE and ACC
# ----------------------------------------------------------------

using NetCDF
using DelimitedFiles
using NCDatasets
using Statistics

function main()
    Panel_1()
    Panel_2()
    Panel_34()
end


function Panel_1()
    fh_ccsm = Dataset(homedir()*"/Data/ECCO/RISE/results/Charleston/Charleston_ccsm_ssh.nc","r")
    fh_obs  = Dataset(homedir()*"/Data/ECCO/RISE/results/Charleston/Charleston_observations.nc","r")
    dir_gmt = homedir()*"/Projects/2020_ECCO_adjoint/GMT/Charleston/supplement/ccsm_drift/"

    t_obs = fh_obs["time"][25:288]
    rsl_obs = fh_obs["TG_lmsl"][25:288]
    rsl_obs = rm_seas(t_obs,rsl_obs)
    t_ccsm = fh_ccsm["tvec"][1,13:276,1] + fh_ccsm["tvec"][1,13:276,2]./12f0 .- 1.f0/24f0
    rsl_ccsm    = mean(fh_ccsm["eta_full"][1,:,13:276],dims=1)[:]
    rsl_dedrift = mean(fh_ccsm["eta_dedrift"][1,:,13:276],dims=1)[:]
    close(fh_ccsm)
    close(fh_obs)

    # Save
    writedlm(dir_gmt*"startpoint.txt",[t_ccsm[1:12:length(t_ccsm)] 100 .* rsl_ccsm[1:12:length(t_ccsm)]],';')
    writedlm(dir_gmt*"noseas_startpoint.txt",[t_ccsm[1:12:length(t_ccsm)] 100 .* rsl_dedrift[1:12:length(t_ccsm)]],';')
    yrs = convert.(Int,unique(floor.(t_obs)))
    for yr in yrs
        acc_idx = floor.(t_obs).==yr
        writedlm(dir_gmt*string(yr)*".txt",[t_ccsm[acc_idx] 100 .* rsl_ccsm[acc_idx]],';')
        writedlm(dir_gmt*string(yr)*"_noseas.txt",[t_ccsm[acc_idx] 100 .* rsl_dedrift[acc_idx]],';')
    end
    writedlm(dir_gmt*"rsl_obs.txt",[t_obs 100 .* rsl_obs],';')
    return nothing
end

function Panel_2()
    dir_gmt = homedir()*"/Projects/2020_ECCO_adjoint/GMT/Charleston/supplement/ccsm_drift/"
    dir_ccsm = homedir()*"/Data/ECCO/RISE/prediction_models/ccsm/raw/JanIC/"

    yrs = [1995:2016...]
    fn_1 = dir_ccsm*"ccsm4_0_cfsrr_Fcst.E1.pop.h.1993-01.nc"
    ϕ=ncread(fn_1,"lon")    
    θ=ncread(fn_1,"lat")
    ydiff = zeros(Float32,length(ϕ),length(θ),length(yrs))
    for ens_member ∈ 1:10
        for (idx,yr) in enumerate(yrs)
            ybegin = ncread(dir_ccsm*"ccsm4_0_cfsrr_Fcst.E"*string(ens_member)*".pop.h."*string(yr)*"-01.nc","SSH")[:,:,1]
            yend   = ncread(dir_ccsm*"ccsm4_0_cfsrr_Fcst.E"*string(ens_member)*".pop.h."*string(yr)*"-12.nc","SSH")[:,:,1]
            @. ybegin[ybegin>1e35] = 0
            @. yend[yend>1e35] = 0
            @. ydiff[:,:,idx] += (yend - ybegin)
        end
    end
    @. ydiff /= 10
    ydiff_mn = mean(ydiff,dims=3)[:,:,1]

    # Write file
    fn_gmt = dir_gmt*"ccsm_drift_map.nc"
    fh = Dataset(fn_gmt,"c")
    defDim(fh,"lon",length(ϕ))
    defDim(fh,"lat",length(θ))
    defVar(fh,"lon",ϕ,("lon",),deflatelevel=4)
    defVar(fh,"lat",θ,("lat",),deflatelevel=4)
    defVar(fh,"ydiff",ydiff_mn,("lon","lat",),deflatelevel=4)
    close(fh)
    return nothing
end

function Panel_34()
    dir_gmt = homedir()*"/Projects/2020_ECCO_adjoint/GMT/Charleston/supplement/ccsm_drift/"
    fh_stats  = Dataset(homedir()*"/Data/ECCO/RISE/results/Charleston/Charleston_ccsm_TG_stats.nc","r")
    writedlm(dir_gmt*"corr_sla.txt",[[1:12...] fh_stats["corr_sla_pred_mean"][:]],';')
    writedlm(dir_gmt*"rmse_sla.txt",[[1:12...] 100 .* fh_stats["rmse_sla_pred_mean"][:]],';')
    writedlm(dir_gmt*"rmse_mss.txt",[[1:12...] 100 .* fh_stats["rmse_obs_seas_mean"][:]],';')
    close(fh_stats)
end

function rm_seas(tval,tseries)
    tseries_noseas = similar(tseries)
    for mnth in 1:12
        mnth_idx = [mnth:12:length(tval)...]
        tseries_noseas[mnth_idx] = tseries[mnth_idx] .- mean(tseries[mnth_idx])
    end
    return tseries_noseas
end
     
