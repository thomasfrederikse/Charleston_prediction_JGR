# ---------------------------------------------------------
# Read all reconstructions and compute RMSE and correlation
# RMSE and correlation for each: 
# - ensemble member
# - forcing component
# - initialziation month
# ---------------------------------------------------------
using LinearAlgebra
using NetCDF
using NCDatasets
using Dates
using Statistics
using DelimitedFiles
using Interpolations
using Distributions
function main()
    settings = def_settings()
    obs = read_obs(settings)
    predictions = read_predictions(obs,settings)
    detrend!(obs,predictions,settings)
    stats = compute_stats(obs,predictions)
    save_data(obs,predictions,stats,settings)

    # Significance of correlation
    # 1. EffDOF
    lag_cor = cor(obs["TG"]["msl_noseas"][2:end],obs["TG"]["msl_noseas"][1:end-1])
    ndof = round(262 * (1-lag_cor)/(1+lag_cor))
    tstat = quantile(Distributions.TDist(ndof),0.95)
    sigcorr = tstat / sqrt(ndof+tstat^2)

    # simple_plot(stats,settings)
end

function def_settings()
    settings=Dict()
    settings["location"] = "Charleston"
    settings["detrend"] = true # Detrend all time series
    settings["obs_truth"] = "TG" # Define which observation is the ground truth for computing RMSE and ACC
    settings["pred_model"] = "ccsm" # Currently: ccsm and spear
    settings["include_IB"] = false
    settings["include_IB"] & (settings["obs_truth"] == "Alt") ? throw("Altimetry and IB: TOO DANGEROUS! DON'T GO THERE!") : nothing
    settings["include_IB"] & (settings["pred_model"] == "ccsm") ? throw("CCSM: no IB ERROR ERROR stop") : nothing
    # Time steps and ensemble members
    settings["months"]   = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
    settings["ens_range"]  = [1:10...] # Array that sets which ensemble members to include
    settings["procs"] = ["empmr","qnet","tauu","tauv"]
    settings["years"]  = [1995:2016...]
    # Observations and prediction methods
    settings["obs"] = ["Alt","TG"] # Observation types
    settings["pred"] = ["obs_seas","sla_pred","sla_obs","conv_obs","conv_seas","conv_obs_pred","conv_obs_seas","persistence"] #The predictors
    settings["pred_ens"] = ["sla_pred","conv_obs_pred"] # Predictors that have ensemble members
    settings["pred_forcing"] = ["conv_obs","conv_seas","conv_obs_seas","conv_obs_pred"] # Predictors that have individual forcings
    # Directories
    settings["dir_data"] = homedir()*"/Data/ECCO/RISE/"
    settings["dir_frc_pred"]     = settings["dir_data"]*"prediction_models/"*settings["pred_model"]*"/regrid/" 
    settings["dir_frc_obs"]     = settings["dir_data"]*"ecco_forcing/"
    settings["dir_adjoint_sens"] = settings["dir_data"]*"adjoint_sensitivities/"*settings["location"]*"/"
    settings["dir_results"] = settings["dir_data"]*"results/"*settings["location"]*"/"
    settings["fn_observations"] = settings["dir_results"]*settings["location"]*"_observations.nc"
    # Coordinates from location file
    settings["ϕ"] = ncread(settings["fn_observations"],"lon")[1]
    settings["θ"] = ncread(settings["fn_observations"],"lat")[1]
    return settings
end

function read_obs(settings)
    # ---------------------------------------------
    # Read all observations (Altimetry, tide gauge)
    # Separate mean seasonal cycle
    # ---------------------------------------------
    println("Reading observations...")
    obs = Dict()
    obs["time"] = [settings["years"][1]+1/24:1/12:settings["years"][end]+1+23/24...]
    obs_time = ncread(settings["fn_observations"],"time")

    # Altimetry
    obs["Alt"] = Dict()
    obs["Alt"]["gmsl"] = LinearInterpolation(obs_time,ncread(settings["fn_observations"],"gmsl"))(obs["time"])
    obs["Alt"]["msl"]  = LinearInterpolation(obs_time,ncread(settings["fn_observations"],"Alt_lmsl"))(obs["time"])

    # Load IB
    obs["TG"] = Dict()
    obs["TG"]["IB"] = LinearInterpolation(obs_time,ncread(settings["fn_observations"],"IB"))(obs["time"])

    if settings["include_IB"]
        obs["IB"] = Dict()
        obs["IB"]["msl"] = obs["TG"]["IB"] 
        obs["TG"]["msl"] = LinearInterpolation(obs_time,ncread(settings["fn_observations"],"TG_lmsl"))(obs["time"]) .+ obs["IB"]["msl"] 
    else
        obs["TG"]["msl"] = LinearInterpolation(obs_time,ncread(settings["fn_observations"],"TG_lmsl"))(obs["time"]) 
    end

    # Separate seasonal cycle
    for obs_idx ∈ settings["obs"]
        obs[obs_idx]["seas"], obs[obs_idx]["msl_seas"], obs[obs_idx]["msl_noseas"] = separate_seas(obs["time"],obs[obs_idx]["msl"])
    end

    if settings["include_IB"]
        obs["IB"]["seas"], obs["IB"]["msl_seas"], obs["IB"]["msl_noseas"] = separate_seas(obs["time"],obs["IB"]["msl"])
    end

    return obs
end

function read_predictions(obs,settings)
    # --------------------------------------------------------------
    # Read predictions
    # - pred SSH
    # - obs SSH / adjoint
    # - Mean seasonal cycle
    # - Hybrid (EECO-CCSM, ECCO-seas)
    #
    # Each projection is corrected as follows:
    # 1. Mean removed
    # 2. Seasonal cycle replaced by:
    #   a. Observed cycle for total sea level
    #   b. Cycle in ECCO reconstruction for individual forcings
    # --------------------------------------------------------------
    println("Reading predictions...")
    predictions = Dict()
    predictions["time"] = [settings["years"][1]+1/24:1/12:settings["years"][end]+1+23/24...]
    predictions["time_stack"] = ts_to_stacked(predictions["time"])
    [predictions[pred] = Dict() for pred ∈ settings["pred"]]

    # Seasonal cycle
    predictions["obs_seas"]["mean"] = ts_to_stacked(obs[settings["obs_truth"]]["msl_seas"])
   
    # Predicted IB effect
    if settings["include_IB"]
        predictions["IB"] = Dict()
        predictions["IB"]["ens"]  = zeros(12,length(settings["ens_range"]),size(predictions["time_stack"],2))
        predictions["IB"]["mean"] = zeros(12,size(predictions["time_stack"],2))
    
        fn = settings["dir_results"] * settings["location"]*"_"*settings["pred_model"]*"_IB.nc"
        pred_IB_time = ncread(fn,"tvec")[:,:,1] + ncread(fn,"tvec")[:,:,2] ./12 .- 1/24
        pred_IB = ncread(fn,"IB_dedrift")
        for mnth ∈ 1:12
            predictions["IB"]["ens"][mnth,:,:] = LinearInterpolation((settings["ens_range"],pred_IB_time[mnth,:]),pred_IB[mnth,:,:],extrapolation_bc=Flat())(settings["ens_range"],predictions["time_stack"][mnth,:])
        end
        predictions["IB"]["mean"][:,:] = mean(predictions["IB"]["ens"],dims=2)[:,1,:]
    end

    # predicted SSH
    predictions["sla_pred"]["ens"]  = zeros(12,length(settings["ens_range"]),size(predictions["time_stack"],2))
    predictions["sla_pred"]["mean"] = zeros(12,size(predictions["time_stack"],2))
    fn = settings["dir_results"] * settings["location"]*"_"*settings["pred_model"]*"_ssh.nc"
    pred_ssh_time = ncread(fn,"tvec")[:,:,1] + ncread(fn,"tvec")[:,:,2] ./12 .- 1/24
    pred_ssh = ncread(fn,"eta_dedrift")
    for mnth ∈ 1:12
        predictions["sla_pred"]["ens"][mnth,:,:] = LinearInterpolation((settings["ens_range"],pred_ssh_time[mnth,:]),pred_ssh[mnth,:,:],extrapolation_bc=Flat())(settings["ens_range"],predictions["time_stack"][mnth,:])
    end
    if settings["include_IB"]
        predictions["sla_pred"]["ens"] += predictions["IB"]["ens"]
    end
    predictions["sla_pred"]["mean"][:,:] = mean(predictions["sla_pred"]["ens"],dims=2)[:,1,:]

    # Persistence
    predictions["persistence"]["mean"] = zeros(12,size(predictions["time_stack"],2))
    fn = settings["dir_results"] * settings["location"]*"_persistence.nc"
    pers_ssh_time = ncread(fn,"tvec")[:,:,1] + ncread(fn,"tvec")[:,:,2] ./12 .- 1/24
    pers_ssh = ncread(fn,"eta_persistence")
    for mnth ∈ 1:12
        predictions["persistence"]["mean"][mnth,:] = LinearInterpolation(pers_ssh_time[mnth,:],pers_ssh[mnth,:], extrapolation_bc=Flat())(predictions["time_stack"][mnth,:])
    end

    # ECCO SLA
    time_ecco, sshdyn_ecco = read_sshdyn(settings)
    predictions["sla_obs"]["mean"] = ts_to_stacked(LinearInterpolation(time_ecco, sshdyn_ecco)(predictions["time"])) 
    if settings["include_IB"]
        predictions["sla_obs"]["mean"] += ts_to_stacked(obs["TG"]["IB"])
    end

    # Hybrid reconstructions and ECCO convolutions
    for pred ∈ settings["pred_forcing"]
        predictions[pred]["proc"] = Dict()
        predictions[pred]["proc"]["mean"] = zeros(length(settings["procs"]),12,size(predictions["time_stack"],2))
        pred in settings["pred_ens"] ? predictions[pred]["proc"]["ens"] = zeros(length(settings["procs"]),12,length(settings["ens_range"]),size(predictions["time_stack"],2)) : nothing
    end

    for mnth ∈ 1:12
        fn = settings["dir_results"]*settings["location"]*"_pred_"*settings["pred_model"]*"_"*settings["months"][mnth]*"IC.nc"
        predictions["conv_obs"]["proc"]["mean"][:,mnth,:] = ncread(fn,"conv_obs")
        predictions["conv_seas"]["proc"]["mean"][:,mnth,:] = ncread(fn,"conv_seas")
        predictions["conv_obs_seas"]["proc"]["mean"][:,mnth,:] = ncread(fn,"conv_obs_seas")
        predictions["conv_obs_pred"]["proc"]["ens"][:,mnth,:,:] = ncread(fn,"conv_obs_pred")
    end

    # Sum over forcing ensemble members to obtain mean forcing
    for pred ∈ settings["pred_forcing"]
        pred in settings["pred_ens"] ? predictions[pred]["proc"]["mean"] = mean(predictions[pred]["proc"]["ens"],dims=3)[:,:,1,:] : nothing
    end 

    # Sum over forcings to obtain sea level
    for pred ∈ settings["pred_forcing"]
        if pred in settings["pred_ens"] 
            predictions[pred]["ens"] = sum(predictions[pred]["proc"]["ens"],dims=1)[1,:,:,:] 
            predictions[pred]["mean"] = mean(predictions[pred]["ens"],dims=2)[:,1,:]
        else
            predictions[pred]["mean"] = sum(predictions[pred]["proc"]["mean"],dims=1)[1,:,:]
        end
    end

    # Add IB
    if settings["include_IB"]
        predictions["conv_obs"]["mean"] += ts_to_stacked(obs["TG"]["IB"])
        predictions["conv_obs_seas"]["mean"] += predictions["IB"]["mean"]
        predictions["conv_obs_pred"]["mean"] += predictions["IB"]["mean"]
        predictions["conv_obs_pred"]["ens"] += predictions["IB"]["ens"]
        predictions["conv_seas"]["mean"] += ts_to_stacked(obs["IB"]["msl_seas"])
    end

    # For all sea-level projections:
    # Remove mean and replace seasonal cycle with observed cycle
    msc_obs = ts_to_stacked(obs[settings["obs_truth"]]["msl_seas"])
    for pred ∈ settings["pred"], mnth ∈ 1:12
        predictions[pred]["mean"][mnth,:] = replace_seas(predictions["time_stack"][mnth,:],predictions[pred]["mean"][mnth,:],msc_obs[mnth,:])
        if haskey(predictions[pred],"ens")
            for ens ∈ settings["ens_range"]
                predictions[pred]["ens"][mnth,ens,:] = replace_seas(predictions["time_stack"][mnth,:],predictions[pred]["ens"][mnth,ens,:],msc_obs[mnth,:])
            end
        end
    end

    # For all projections with individual forcings: remove mean and replace with ECCO convolution forcing seasonal cycle
    for pred ∈ settings["pred_forcing"], mnth ∈ 1:12, proc ∈ 1:4
        predictions[pred]["proc"]["mean"][proc,mnth,:] = replace_seas(predictions["time_stack"][mnth,:],predictions[pred]["proc"]["mean"][proc,mnth,:],predictions["conv_seas"]["proc"]["mean"][proc,mnth,:])
        if haskey(predictions[pred]["proc"],"ens")
            for ens ∈ settings["ens_range"]
                predictions[pred]["proc"]["ens"][proc,mnth,ens,:] = replace_seas(predictions["time_stack"][mnth,:],predictions[pred]["proc"]["ens"][proc,mnth,ens,:],predictions["conv_seas"]["proc"]["mean"][proc,mnth,:])
            end
        end
    end
    return predictions
end

function detrend!(obs,predictions,settings)
    if settings["detrend"]
        for obs_idx ∈ settings["obs"]
            obs[obs_idx]["msl"] = detrend(obs["time"],obs[obs_idx]["msl"])
        end
        for pred ∈ settings["pred"] 
            for mnth ∈ 1:12
                predictions[pred]["mean"][mnth,:] = detrend(predictions["time_stack"][mnth,:],predictions[pred]["mean"][mnth,:])
                haskey(predictions[pred],"ens") ? [predictions[pred]["ens"][mnth,ens,:] = detrend(predictions["time_stack"][mnth,:],predictions[pred]["ens"][mnth,ens,:]) for ens ∈ settings["ens_range"]] : nothing
            end
        end
    end
    return nothing
end

function compute_stats(obs,predictions)
    # ---------------------------------------------------------------------------
    # Compute statistics from time series
    # For each prediction (ensemble mean and individual ensemble members):
    #  - RMSE with respect to observations
    #  - ACC with respect to observations and seasonal cycle
    # For each prediction with individual forcings:
    #  - RMSE with respect to ECCO convolution
    # ---------------------------------------------------------------------------
    println("Computing statistics...")
    obs_stacked      = ts_to_stacked(obs[settings["obs_truth"]]["msl"])
    obs_seas_stacked = ts_to_stacked(obs[settings["obs_truth"]]["msl_seas"])
    stats = Dict()
    # Total sea level
    for pred ∈ settings["pred"]
        stats[pred] = Dict()
        stats[pred]["rmse_mean"] = compute_rmse(obs_stacked,predictions[pred]["mean"]) 
        stats[pred]["corr_mean"] = compute_corr(obs_stacked,predictions[pred]["mean"],obs_seas_stacked)
        if haskey(predictions[pred],"ens")
            stats[pred]["rmse_ens"] = zeros(12,length(settings["ens_range"]))
            stats[pred]["corr_ens"] = zeros(12,length(settings["ens_range"]))
            for ens ∈ settings["ens_range"]
                stats[pred]["rmse_ens"][:,ens] = compute_rmse(obs_stacked,predictions[pred]["ens"][:,ens,:]) 
                stats[pred]["corr_ens"][:,ens] = compute_corr(obs_stacked,predictions[pred]["ens"][:,ens,:],obs_seas_stacked)
            end
        end
    end

    # Forcings
    for pred ∈ settings["pred_forcing"]
        stats[pred]["proc"] = Dict()
        stats[pred]["proc"]["rmse_mean"] = zeros(4,12)
        if haskey(predictions[pred]["proc"],"ens")
            stats[pred]["proc"]["rmse_ens"] = zeros(4, 12,length(settings["ens_range"]))
        end
        for proc ∈ 1:4
            stats[pred]["proc"]["rmse_mean"][proc,:] = compute_rmse(predictions[pred]["proc"]["mean"][proc,:,:],predictions["conv_obs"]["proc"]["mean"][proc,:,:])
            if haskey(predictions[pred]["proc"],"ens")
                for ens ∈ settings["ens_range"]
                    stats[pred]["proc"]["rmse_ens"][proc,:,ens] = compute_rmse(predictions[pred]["proc"]["ens"][proc,:,ens,:],predictions["conv_obs"]["proc"]["mean"][proc,:,:])
                end
            end
        end
    end

    # Corr and RMSE as function of initialization month
    for pred ∈ settings["pred"]
        stats[pred]["rmse_mean_init_seas"] = compute_rmse_init_seas(obs_stacked,predictions[pred]["mean"]) 
        stats[pred]["corr_mean_init_seas"] = compute_acc_init_seas(obs_stacked,predictions[pred]["mean"],obs_seas_stacked)
    end

    return stats
end

# -----------------
# Generic functions
# -----------------
function detrend(time,tseries)
    amat = ones(length(time),6);        
    amat[:,3] = sin.(2*pi.*time)
    amat[:,4] = sin.(4*pi.*time)
    amat[:,5] = cos.(2*pi.*time)
    amat[:,6] = cos.(4*pi.*time)    
    amat[:,2] = time.- mean(time);
    sol = amat\tseries
    sol[3:end] .= 0
    tseries_detrend = tseries - amat * sol
    return tseries_detrend
end

function separate_seas(time,tseries)
    tseries_detrend = detrend(time,tseries)
    seas = zeros(12)
    tseries_seas = similar(tseries)
    for mnth ∈ 1:12
        mnth_idx = [mnth:12:length(time)...]
        seas[mnth] = mean(tseries_detrend[mnth_idx])
        tseries_seas[mnth_idx] .= seas[mnth]
    end
    tseries_noseas = tseries - tseries_seas
    return seas, tseries_seas, tseries_noseas
end

function ts_to_stacked(tseries)
    # From a time series of length n , make a stacked time series, where the first row runs from 1:n-12, second 2:n-11 etc
    tseries_stacked = zeros(12,length(tseries)-12)
    for mnth ∈ 1:12
        @. tseries_stacked[mnth,:] = tseries[mnth:end-13+mnth]
    end
    return tseries_stacked
end

function replace_seas(time,tseries,seas)
    tseries_detrend = detrend(time,tseries)
    tseries_seas = similar(tseries)
    for mnth ∈ 1:12
        mnth_idx = [mnth:12:length(time)...]
        tseries_seas[mnth_idx] .= mean(tseries_detrend[mnth_idx])
    end
    tseries = tseries - tseries_seas + seas
    tseries .-= mean(tseries)
    return tseries
end

function compute_rmse_init_seas(truth,model)
    # rmse_init_mnth = zeros(Float32,12,12)
    # for init_mnth in 1:12
    #     for lag_mnth in 1:12
    #         rmse_init_mnth[init_mnth,lag_mnth] = sqrt(mean((model[init_mnth,lag_mnth:12:end] - truth[init_mnth,lag_mnth:12:end]).^2))
    #     end
    # end
    init_seas_t = zeros(Int,4,3)
    init_seas_t[1,:] = [12,1,2]
    init_seas_t[2,:] = [3,4,5]
    init_seas_t[3,:] = [6,7,8]
    init_seas_t[4,:] = [9,10,11]
    rmse_init_seas = zeros(Float32,4,12)
    for init_seas in 1:4
        for lag_mnth in 1:12
            rmse_init_seas[init_seas,lag_mnth] = sqrt(mean((model[init_seas_t[init_seas,:],lag_mnth:12:end] - truth[init_seas_t[init_seas,:],lag_mnth:12:end]).^2))
        end
    end
    return rmse_init_seas
end

function compute_acc_init_seas(truth,model,seasonal)
    init_seas_t = zeros(Int,4,3)
    init_seas_t[1,:] = [12,1,2]
    init_seas_t[2,:] = [3,4,5]
    init_seas_t[3,:] = [6,7,8]
    init_seas_t[4,:] = [9,10,11]

    acc_init_seas = zeros(Float32,4,12)
    for init_seas in 1:4
        for lag_mnth in 1:12
            acc_init_seas[init_seas,lag_mnth] = cor(model[init_seas_t[init_seas,:],lag_mnth:12:end][:]-seasonal[init_seas_t[init_seas,:],lag_mnth:12:end][:],truth[init_seas_t[init_seas,:],lag_mnth:12:end][:]-seasonal[init_seas_t[init_seas,:],lag_mnth:12:end][:])
        end
    end
    return acc_init_seas
end


function compute_rmse(truth,model)
    rmse = zeros(12)
    for mnth=1:12
        mnth_idx = [mnth:12:264...]
        model_lcl = model[:,mnth_idx][:]
        truth_lcl = truth[:,mnth_idx][:]
        rmse[mnth] = sqrt(mean((model_lcl - truth_lcl).^2))
    end
    return rmse
end

function compute_corr(truth,model,seasonal)
    corr = zeros(12)
    for mnth=1:12
        mnth_idx = [mnth:12:264...]
        model_lcl = model[:,mnth_idx][:]
        truth_lcl = truth[:,mnth_idx][:]
        seasonal_lcl = seasonal[:,mnth_idx][:]
        corr[mnth] = cor(model_lcl-seasonal_lcl,truth_lcl-seasonal_lcl)
    end
    return corr
end

function read_sshdyn(settings)
    # Read ECCO info
    dir_ecco_sshdyn = homedir() * "/Data/ECCO/v4r4/SSHDYN/"
    ϕ_ECCO = ncread(dir_ecco_sshdyn * "1992/SSHDYN_1992_01.nc","XC")
    θ_ECCO = ncread(dir_ecco_sshdyn * "1992/SSHDYN_1992_01.nc","YC")
    @. ϕ_ECCO[ϕ_ECCO.<0] += 360
    coords = argmin(@. ((settings["ϕ"]-ϕ_ECCO)^2 + (settings["θ"]-θ_ECCO)^2))
    area = ncread(homedir() * "/Data/ECCO/v4r4/ECCOv4r4_grid.nc","rA")
    slm = ncread(homedir() * "/Data/ECCO/v4r4/ECCOv4r4_grid.nc","hFacW",start=[1,1,1,1],count=[-1,-1,-1,1])[:,:,:,1];
    area_tot = sum(area.*slm)

    time_ecco = [1992+1/24:1/12:2017+23/24...]
    tvec_ecco = zeros(Int,length(time_ecco),2)
    @. tvec_ecco[:,1] = floor(time_ecco)
    @. tvec_ecco[:,2] = round((time_ecco - tvec_ecco[:,1]) * 12 + 0.5)

    ssh_lcl = zeros(Float32,length(time_ecco))
    ssh_glb = zeros(Float32,length(time_ecco))

    # Read data
    for tstep ∈ 1:length(time_ecco)
        fn = dir_ecco_sshdyn * string(tvec_ecco[tstep,1])*"/SSHDYN_"*string(tvec_ecco[tstep,1])*"_"*lpad(tvec_ecco[tstep,2],2,"0")*".nc"
        ssh_ts = dropdims(ncread(fn,"SSHDYN"),dims=4);
        ssh_lcl[tstep] = ssh_ts[coords]
        ssh_glb[tstep] = sum(area.*ssh_ts) / area_tot
    end

    # Save to standard format
    sshdyn_ecco = ssh_lcl - ssh_glb
    return time_ecco, sshdyn_ecco
end

function save_data(obs,predictions,stats,settings)
    # Save time series and statistics to 2 netCDF files
    println("Saving...")
    settings["detrend"] ? fn_ts=settings["dir_results"]*settings["location"]*"_"*settings["pred_model"]*"_"*settings["obs_truth"]*"_detrend_tseries.nc" : fn_ts=settings["dir_results"]*settings["location"]*"_"*settings["pred_model"]*"_"*settings["obs_truth"]*"_tseries.nc"
    settings["detrend"] ? fn_stats=settings["dir_results"]*settings["location"]*"_"*settings["pred_model"]*"_"*settings["obs_truth"]*"_detrend_stats.nc" : fn_stats=settings["dir_results"]*settings["location"]*"_"*settings["pred_model"]*"_"*settings["obs_truth"]*"_stats.nc"
    
    # DateTime
    time_array = @. DateTime(convert(Int,floor(predictions["time_stack"])), round((predictions["time_stack"] - floor(predictions["time_stack"]))*12+0.5),15)
    tstart_array = repeat(predictions["time_stack"][:,1:12:end],inner=[1,12])
    tstart_array = @. DateTime(convert(Int,floor(tstart_array)), round((tstart_array - floor(tstart_array))*12+0.5),15)

    # Save time series into netCDF format
    fh = Dataset(fn_ts,"c")
    defDim(fh,"tsteps",size(predictions["time_stack"],2))
    defDim(fh,"tstart",size(predictions["time_stack"],1))
    defDim(fh,"proc",length(settings["procs"]))
    defDim(fh,"ens",length(settings["ens_range"]))

    defVar(fh,"time",time_array,("tstart","tsteps"),deflatelevel=4)
    defVar(fh,"tstart",tstart_array,("tstart","tsteps"),deflatelevel=4)
    defVar(fh,"ens",settings["ens_range"],("ens",),deflatelevel=4)
    defVar(fh,"proc",settings["procs"],("proc",),deflatelevel=4)
    
    for obs_idx ∈ settings["obs"]
        defVar(fh,"obs_"*obs_idx,ts_to_stacked(obs[obs_idx]["msl"]),("tstart","tsteps"),deflatelevel=4,attrib = Dict("units" => "m","long_name" => "Observed mean sea level ("*obs_idx*")",))
        defVar(fh,"obs_"*obs_idx*"_seas",ts_to_stacked(obs[obs_idx]["msl_seas"]),("tstart","tsteps"),deflatelevel=4,attrib = Dict("units" => "m","long_name" => "Observed seasonal cycle ("*obs_idx*")",))
    end

    for pred ∈ settings["pred"]
        defVar(fh,"pred_"*pred*"_mean",predictions[pred]["mean"],("tstart","tsteps",),deflatelevel=4,attrib = Dict("units" => "m","long_name" => "Sea-level prediction "*pred))
        haskey(predictions[pred],"ens") ? defVar(fh,"pred_"*pred*"_ens",predictions[pred]["ens"],("tstart","ens","tsteps",),deflatelevel=4,attrib = Dict("units" => "m","long_name" => "Sea-level prediction "*pred* "ensemble members")) : nothing
        if haskey(predictions[pred],"proc")
            defVar(fh,"pred_"*pred*"_proc_mean",predictions[pred]["proc"]["mean"],("proc","tstart","tsteps",),deflatelevel=4,attrib = Dict("units" => "m","long_name" => "Sea-level prediction "*pred*": Individual processes"))
            haskey(predictions[pred]["proc"],"ens") ? defVar(fh,"pred_"*pred*"_proc_ens",predictions[pred]["proc"]["ens"],("proc","tstart","ens","tsteps",),deflatelevel=4,attrib = Dict("units" => "m","long_name" => "Sea-level prediction "*pred* ": Individual processes ensemble members")) : nothing
        end
    end
    close(fh)

    # Save statistics into netCDF format
    fh = Dataset(fn_stats,"c")
    defDim(fh,"months",12)
    defDim(fh,"proc",length(settings["procs"]))
    defDim(fh,"ens",length(settings["ens_range"]))
    defDim(fh,"seasons",4)

    defVar(fh,"months",[1:12...],("months",),deflatelevel=4)
    defVar(fh,"ens",settings["ens_range"],("ens",),deflatelevel=4)
    defVar(fh,"proc",settings["procs"],("proc",),deflatelevel=4)
    defVar(fh,"seasons",["DJF", "MAM", "JJA", "SON"],("seasons",),deflatelevel=4)
 
    for pred ∈ settings["pred"]
        defVar(fh,"rmse_"*pred*"_mean",stats[pred]["rmse_mean"],("months",),deflatelevel=4,attrib = Dict("units" => "m","long_name" => "Mean RMSE of "*pred*" with respect to observed sea level ("*settings["obs_truth"]*")"))
        defVar(fh,"corr_"*pred*"_mean",stats[pred]["corr_mean"],("months",),deflatelevel=4,attrib = Dict("units" => "m","long_name" => "Mean ACC of "*pred*" with respect to observed sea level ("*settings["obs_truth"]*")"))
        if haskey(stats[pred],"rmse_ens")
            defVar(fh,"rmse_"*pred*"_ens",stats[pred]["rmse_ens"],("months","ens",),deflatelevel=4,attrib = Dict("units" => "m","long_name" => "Ensemble RMSE of "*pred*" with respect to observed sea level ("*settings["obs_truth"]*")"))
            defVar(fh,"corr_"*pred*"_ens",stats[pred]["corr_ens"],("months","ens",),deflatelevel=4,attrib = Dict("units" => "m","long_name" => "Ensemble ACC of "*pred*" with respect to observed sea level ("*settings["obs_truth"]*")"))
        end
        if haskey(stats[pred],"proc")
            defVar(fh,"rmse_"*pred*"_proc_mean",stats[pred]["proc"]["rmse_mean"],("proc","months",),deflatelevel=4,attrib = Dict("units" => "m","long_name" => "Individual process mean RMSE of "*pred*" with respect ECCO convolution"))
            haskey(stats[pred]["proc"],"rmse_ens") ? defVar(fh,"rmse_"*pred*"_proc_ens",stats[pred]["proc"]["rmse_ens"],("proc","months","ens"),deflatelevel=4,attrib = Dict("units" => "m","long_name" => "Individual process ensemble RMSE of "*pred*" with respect ECCO convolution")) : nothing
        end
        # Seasonal averages
        defVar(fh,"rmse_"*pred*"_mean_init_seas",stats[pred]["rmse_mean_init_seas"],("seasons","months",),deflatelevel=4,attrib = Dict("units" => "m","long_name" => "Mean RMSE of "*pred*" with respect to observed sea level ("*settings["obs_truth"]*") for each initialization season"))
        defVar(fh,"corr_"*pred*"_mean_init_seas",stats[pred]["corr_mean_init_seas"],("seasons","months",),deflatelevel=4,attrib = Dict("units" => "-","long_name" => "Mean ACC of "*pred*" with respect to observed sea level ("*settings["obs_truth"]*") for each initialization season"))
    end


    close(fh)
    return nothing
end

function simple_plot(stats,settings)
    # RMSE Plots
    plot(stats["obs_seas"]["rmse_mean"],label="Seasonal cycle",linewidth=2)
    plot!(stats["sla_obs"]["rmse_mean"],label="ECCO sla",linewidth=2)
    plot!(stats["conv_obs_pred"]["rmse_mean"],label="Hybrid: ECCO + "*settings["pred_model"],linewidth=2)
    plot!(stats["conv_obs_seas"]["rmse_mean"],label="Hybrid: ECCO + seasonal",linewidth=2)
    plot!(stats["conv_obs"]["rmse_mean"],label="ECCO convolution",linewidth=2)
    plot!(stats["sla_pred"]["rmse_mean"],label=settings["pred_model"]*" sla",linewidth=2,xticks=(1:3:12))
    ylabel!("RMSE (m)")
    xlabel!("Months after initialization",legend=:bottomright)
    title!(settings["location"]*"  "*settings["pred_model"])
    savefig(settings["location"]*"_"*settings["pred_model"]*"_IB_rmse_2yr.png")

    # Correlation Plots
    plot(stats["obs_seas"]["corr_mean"],label="Seasonal cycle",linewidth=2)
    plot!(stats["sla_obs"]["corr_mean"],label="ECCO sla",linewidth=2)
    plot!(stats["conv_obs_pred"]["corr_mean"],label="Hybrid: ECCO + "*settings["pred_model"],linewidth=2)
    plot!(stats["conv_obs_seas"]["corr_mean"],label="Hybrid: ECCO + seasonal",linewidth=2)
    plot!(stats["conv_obs"]["corr_mean"],label="ECCO convolution",linewidth=2)
    plot!(stats["sla_pred"]["corr_mean"],label=settings["pred_model"]*" sla",legend=:bottomright,linewidth=2,xticks=(1:3:12))
    ylabel!("Correlation (m)")
    xlabel!("Months after initialization",legend=:topright)
    title!(settings["location"]*"  "*settings["pred_model"])
    savefig(settings["location"]*"_"*settings["pred_model"]*"_IB_corr_2yr.png")

    # Test initialization month
    cpal = palette(:tab10)
    plot(mean(stats["conv_obs_seas"]["rmse_mean_init_month"][[6,7,8],1:12],dims=1)[:],color=cpal[1])
    plot!(mean(stats["conv_obs_pred"]["rmse_mean_init_month"][[6,7,8],1:12],dims=1)[:],ls=:dash,color=cpal[1])

    plot!(mean(stats["conv_obs_seas"]["rmse_mean_init_month"][[9,10,11],1:12],dims=1)[:],color=cpal[2])
    plot!(mean(stats["conv_obs_pred"]["rmse_mean_init_month"][[9,10,11],1:12],dims=1)[:],ls=:dash,color=cpal[2])

    plot!(mean(stats["conv_obs_seas"]["rmse_mean_init_month"][[12,1,2],1:12],dims=1)[:],color=cpal[3])
    plot!(mean(stats["conv_obs_pred"]["rmse_mean_init_month"][[12,1,2],1:12],dims=1)[:],ls=:dash,color=cpal[3])

    plot!(mean(stats["conv_obs_seas"]["rmse_mean_init_month"][[3,4,5],1:12],dims=1)[:],color=cpal[4])
    plot!(mean(stats["conv_obs_pred"]["rmse_mean_init_month"][[3,4,5],1:12],dims=1)[:],ls=:dash,color=cpal[4])
end
    

function IB_predictability()
    corr_IB = compute_corr(ts_to_stacked(obs["IB"]["msl"]),predictions["IB"]["mean"],ts_to_stacked(obs["IB"]["msl_seas"]))
    rmse_spear = compute_rmse(ts_to_stacked(obs["IB"]["msl"]),predictions["IB"]["mean"])
    rmse_mss   = compute_rmse(ts_to_stacked(obs["IB"]["msl"]),ts_to_stacked(obs["IB"]["msl_seas"]))

    plot(corr_IB,label=settings["pred_model"],legend=:topright,linewidth=2,xticks=(1:3:12))
    title!("Correlation with ERA5 IB")
    savefig("corr_IB.png")

    plot(rmse_spear,label=settings["pred_model"],legend=:topright,linewidth=2,xticks=(1:3:12))
    plot!(rmse_mss,label="Seasonal cycle",legend=:topleft,linewidth=2,xticks=(1:3:12))
    title!("RMSE w.r.t. ERA5 IB")
    savefig("rmse_IB.png")

    plot(obs["time"],obs["IB"]["msl_noseas"],label="IB (ERA5)",linewidth=2)
    plot!(obs["time"],obs["TG"]["msl_noseas"],label="Observations (TG)",linewidth=2)
    ylabel!("Sea level (m)")
    savefig("tg_noseas.png")

    plot(obs["time"],obs["IB"]["msl"],label="IB (ERA5)",linewidth=2)
    plot!(obs["time"],obs["TG"]["msl"],label="Observations (TG)",linewidth=2)
    ylabel!("Sea level (m)")
    savefig("tg_seas.png")
end
