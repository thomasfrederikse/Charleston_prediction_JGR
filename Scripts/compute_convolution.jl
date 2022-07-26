# --------------------------------------------------------------
# Reconstruct sea level using adjoint sensitivities from ECCO
#
# Three reconstructions:
# - "obs"      : convoluition with ECCO/observed forcings
# - "obs_pred" : convolution with observed and predicted forcing
# - "obs_seas" : convolution with observed and seasonal forcing
#
# Forcing from
# obs: ECCO (observation based)
# pred: ensemble forecasts 
# --------------------------------------------------------------
using LinearAlgebra
using NCDatasets
using LoopVectorization
using Dates
using Statistics
using NetCDF
function main()
    settings = def_settings()
    for ic_mnth in settings["month_start"]
        println("IC month "*string(ic_mnth)*"...")
        conv_obs,conv_seas,conv_obs_seas,conv_obs_pred = convolve_all(settings,ic_mnth);
        save_reconstruction(conv_obs,conv_seas,conv_obs_seas,conv_obs_pred,ic_mnth,settings);
    end
end

function def_settings()
    println("Settings...")
    settings = Dict()
    settings["location"] = "Charleston"
    settings["pred_model"] = "ccsm" # Currently: ccsm and spear
    settings["exclude_pred_mnth_in_bc"] = true # To make reviewer 2 happy switch to true
    settings["months"]   = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
    settings["month_start"] = [1:1...] # Start months to include (12 total)
    settings["ens_range"]  = [1:1...] # Array that sets which ensemble members to include
    settings["processes"] = ["empmr","qnet","tauu","tauv"] # Set the processes to be included in the convolution
    settings["processes_long"] = ["Freshwater flux","Heat flux","Zonal wind","Meridional wind"] 
    # Times and dates
    settings["years"] = [1995:2016...] # Years for which the forecast is made. Counted from the initialization date (i.e. initialization in June 1993 provides [June 1993-May 1994] forecast)
    settings["time"] = Array{Date}(undef,length(settings["month_start"]),12*length(settings["years"])) # Array with time steps for each initialization month
    settings["time_start"] = similar(settings["time"]) # Initialization time for each time step for each initialization month
    for mnth in 1:length(settings["month_start"])
        settings["time"][mnth,:] = [Date(settings["years"][1],settings["month_start"][mnth],15):Month(1):Date(2022,1,15)...][1:(12*length(settings["years"]))]
        settings["time_start"][mnth,:] = repeat([Date(yr,settings["month_start"][mnth],1) for yr in settings["years"]],inner=12)
    end
    settings["month"] = [Month(t).value for t in settings["time"]] # Number of month of each time step
return settings
end

function convolve_all(settings,ic_mnth)
    println("Convolving...")
    # ---------------------------------------------------------------
    # Do the actual convolution of the sensitivities with the forcing
    # ---------------------------------------------------------------

    # Prediction forcing time  
    fn_init  = homedir()*"/Data/ECCO/RISE/prediction_models/"*settings["pred_model"]*"/regrid/"*settings["pred_model"]*"_init_"*lpad(ic_mnth,2,"0")*"_ens_01.nc"
    pred_tvec = convert.(Int32,ncread(fn_init,"tvec"))
    pred_t    = @. Date(pred_tvec[:,1],pred_tvec[:,2],15);
    pred_acc  = findall(@. (pred_t-settings["time"][ic_mnth,1]>-Day(800)) & (pred_t- settings["time"][ic_mnth,end] < Day(50)))
    pred_tvec = pred_tvec[pred_acc,:]
    pred_t    = pred_t[pred_acc]

    # Observations forcing time 
    obs_t  = zeros(Date,1356)
    obs_hr = zeros(Int32,1356)
    [obs_t[t] = t * Day(7) + Date(1992,1,1) for t ∈ 1:1356]
    [obs_hr[t] = (Hour(obs_t[t]-Date(1992,1,1))).value for t ∈ 1:1356]
    @. obs_t .-= Day(3)
    obs_acc  = findall(@. (obs_t-settings["time"][ic_mnth,1]>-Day(800)) & (obs_t- settings["time"][ic_mnth,end] < Day(50)))
    obs_t  = obs_t[obs_acc]
    obs_hr = obs_hr[obs_acc]
    
    # pre-allocate forcing and sensitivity
    obs_frc       = Array{Float32}(undef,90,90,13,length(obs_t));
    seas_frc      = Array{Float32}(undef,90,90,13,length(obs_t));
    pred_frc      = Array{Float32}(undef,90,90,13,length(settings["ens_range"]),length(pred_t));
    pred_frc_corr = Array{Float32}(undef,90,90,13,length(settings["ens_range"]),length(pred_t));
    sensitivity = Array{Sensitivity,1}(undef,12)

    # Pre-allocate convolution arrays
    conv_obs      = zeros(Float32,length(settings["processes"]),size(settings["time"],2));  # Convolution with ECCO forcing
    conv_seas     = zeros(Float32,length(settings["processes"]),size(settings["time"],2));  # Convolution with the mean seasonal cycle of the ECCO forcing
    conv_obs_seas = zeros(Float32,length(settings["processes"]),size(settings["time"],2));  # Convolution with mean seasonal cycle of the ECCO forcing for predictive adjoint time steps and ECCO forcing for observed adjoint time steps
    conv_obs_pred = zeros(Float32,length(settings["processes"]),length(settings["ens_range"]),size(settings["time"],2)); # Convolution with bias-corrected CCSM forcing for predictive adjoint time steps and ECCO forcing for observed adjoint time steps
    
    # Loop over each forcing
    for proc_idx in 1:length(settings["processes"]) 
        println("  Convolving "*settings["processes_long"][proc_idx]*"...")
        proc_name = settings["processes"][proc_idx]

        read_pred_forcing!(proc_name,ic_mnth,pred_acc,pred_frc,settings); # Read predict forcing ensemble
        read_ecco_forcing!(proc_name,obs_hr,obs_frc);                     # Read ECCO forcing
        read_sensitivity!(sensitivity,proc_name,settings);                # Read adjoint sensitivity

        if settings["exclude_pred_mnth_in_bc"] == false
            comp_seas_forcing!(obs_frc,seas_frc,obs_t);                       # Compute mean seasonal cycle of ECCO forcing
            corr_forcing_bias!(pred_frc_corr,obs_t,obs_frc,pred_t,pred_frc);                # Apply bias correction to CCSM forcing
        end

        println("    Convolution...")
        Threads.@threads for tstep_pred in 1:size(settings["time"],2)  # Loop over each time step where sea level must be predicted [tstep_pred]
            println(tstep_pred)
            # Determine time steps
            time_sens = @. (sensitivity[settings["month"][ic_mnth,tstep_pred]].time  - Year(2004) + Year(settings["time"][ic_mnth,tstep_pred])) # Determine sensitivity time steps to be used
            time_pred = @. (time_sens > settings["time_start"][ic_mnth,tstep_pred]) # Time steps of the adjoint sensitivity where predicted forcing must be used (other time step use observed forcing)
            # Pre-allocate some arrays
            adj_conv_pred    = zeros(Float32,length(settings["ens_range"]))
            sens_lcl         = zeros(Float32,90,90,13);

            ### Bias correct using R2 approach
            if settings["exclude_pred_mnth_in_bc"] == true
                # Determine obs and pred forcings to exclude
                tsteps_obs_inc = ones(Bool,length(obs_t))
                tsteps_pred_inc = ones(Bool,length(pred_t))
                for exc_idx in 1:length(time_sens) # Loop over each sensitivity time step [tstep_sens] belonging to [tstep_pred]
                    tsteps_obs_inc[argmin(abs.(obs_t - time_sens[exc_idx]))] = false # Determine the ECCO forcing time step corresponding to the sensitivity time step
                    tsteps_pred_inc[findfirst(pred_t .== Date(Year(time_sens[exc_idx]).value,Month(time_sens[exc_idx]).value,15))] = false; # Determine the prediction forcing time step corresponding to the sensitivity time step
                end
                # Seasonal forcings
                comp_seas_forcing_no_ovl!(obs_frc,seas_frc,obs_t,tsteps_obs_inc)
                corr_forcing_bias_no_ovl!(pred_frc_corr,obs_t,obs_frc,pred_t,pred_frc,tsteps_obs_inc)
            end

            for idx_sens in 1:length(time_sens) # Loop over each sensitivity time step [tstep_sens] belonging to [tstep_pred]
                @. sens_lcl = sensitivity[settings["month"][ic_mnth,tstep_pred]].sens[:,:,:,idx_sens];
                tstep_obs_frc  = argmin(abs.(obs_t - time_sens[idx_sens])) # Determine the ECCO forcing time step corresponding to the sensitivity time step
                tstep_pred_frc = findfirst(pred_t .== Date(Year(time_sens[idx_sens]).value,Month(time_sens[idx_sens]).value,15)); # Determine the prediction forcing time step corresponding to the sensitivity time step
                
                # Sea-level contribution for [tstep_pred] from convolving forcing and sensitivity at [tstep_sens] 
                adj_conv_obs  = vmapreduce(*,+,sens_lcl,@views obs_frc[:,:,:,tstep_obs_frc])
                adj_conv_seas = vmapreduce(*,+,sens_lcl,@views seas_frc[:,:,:,tstep_obs_frc])
                for ens in 1:length(settings["ens_range"])
                    adj_conv_pred[ens] = vmapreduce(*,+,sens_lcl,@views(pred_frc_corr[:,:,:,ens,tstep_pred_frc]))
                end

                conv_obs[proc_idx,tstep_pred]  += adj_conv_obs
                conv_seas[proc_idx,tstep_pred] += adj_conv_seas
                conv_obs_seas[proc_idx,tstep_pred]      += ifelse(time_pred[idx_sens],adj_conv_seas,adj_conv_obs)
                @. conv_obs_pred[proc_idx,:,tstep_pred] += ifelse(time_pred[idx_sens],adj_conv_pred,adj_conv_obs)
            end
        end
    end
    return conv_obs,conv_seas,conv_obs_seas,conv_obs_pred
end

function read_pred_forcing!(proc_name,ic_mnth,pred_acc,pred_frc,settings)
    println("    Reading prediction forcing...")
    # Forcing time
    for ens in 1:length(settings["ens_range"])
        fn = homedir()*"/Data/ECCO/RISE/prediction_models/"*settings["pred_model"]*"/regrid/"*settings["pred_model"]*"_init_"*lpad(ic_mnth,2,"0")*"_ens_"*lpad(ens,2,"0")*".nc"
        pred_frc[:,:,:,ens,:] = ncread(fn,proc_name,start=[1,1,1,pred_acc[1]],count=[-1,-1,-1,length(pred_acc)]);
    end
    return nothing
end

function read_ecco_forcing!(proc_name,obs_hr,obs_frc)
    println("    Reading ECCO forcing...")
    dir_forcing = homedir()*"/Data/ECCO/RISE/ecco_forcing/"
    Threads.@threads for idx ∈ 1:length(obs_hr)
        filename = dir_forcing * proc_name * "_weekly_v1."*lpad(string(obs_hr[idx]),10,'0')*".data"
        read_forcing_file!(obs_frc, idx, filename)
    end
    return nothing
end

function comp_seas_forcing!(obs_frc,seas_frc,obs_t)
    # Compute the seasonal cycle of the forcing
    tval = [(obs_t[i]-obs_t[1]).value for i ∈ 1:length(obs_t)]./365.24
    amat = ones(Float32,length(tval),10)
    amat[:,2] = tval.-mean(tval)
    @. amat[:,3] = sin(2*π*tval)
    @. amat[:,4] = cos(2*π*tval)
    @. amat[:,5] = sin(4*π*tval)
    @. amat[:,6] = cos(4*π*tval)
    @. amat[:,7] = sin(6*π*tval)
    @. amat[:,8] = cos(6*π*tval)
    @. amat[:,9] = sin(8*π*tval)
    @. amat[:,10] = cos(8*π*tval)
    amat_tr = transpose(amat);
    amat_sq = inv(amat_tr*amat);
    Threads.@threads for i ∈ 1:size(obs_frc,1)
        am_m_ts = zeros(Float32,10)
        sol = zeros(Float32,10)
        ts_sol = zeros(Float32,length(tval))
        @inbounds for j ∈ 1:size(obs_frc,2)
            @simd for k ∈ 1:size(obs_frc,3)
                mul!(am_m_ts,amat_tr, @views obs_frc[i,j,k,:]);
                mul!(sol,amat_sq,am_m_ts);
                sol[2] = 0 # Do not reconstruct trend (that'd be cheating ;)
                mul!(ts_sol,amat,sol)
                @. seas_frc[i,j,k,:] = ts_sol;
            end
        end
    end
    return nothing
end

function comp_seas_forcing_no_ovl!(obs_frc,seas_frc,obs_t,tsteps_obs_inc)
    # Compute the seasonal cycle of the forcing
    tval = [(obs_t[i]-obs_t[1]).value for i ∈ 1:length(obs_t)]./365.24

    # Amat for estimation
    amat = ones(Float32,length(tval),10)
    @. amat[:,2] = tval - $mean(tval)
    @. amat[:,3] = sin(2*π*tval)
    @. amat[:,4] = cos(2*π*tval)
    @. amat[:,5] = sin(4*π*tval)
    @. amat[:,6] = cos(4*π*tval)
    @. amat[:,7] = sin(6*π*tval)
    @. amat[:,8] = cos(6*π*tval)
    @. amat[:,9] = sin(8*π*tval)
    @. amat[:,10] = cos(8*π*tval)

    amat_tr = copy(transpose(amat[tsteps_obs_inc,:]));
    amat_sq = inv(amat_tr*amat[tsteps_obs_inc, :])
    for k in 1:size(obs_frc,3)
        @inbounds for j ∈ 1:size(obs_frc,2)
            for i ∈ 1:size(obs_frc,1)          
                sol = amat_sq * (amat_tr * obs_frc[i,j,k,tsteps_obs_inc])
                sol[2] = 0.0f0 # Do not reconstruct trend (that'd be cheating ;)
                ts_sol = amat * sol
                @turbo for t in 1:length(obs_t)
                    seas_frc[i,j,k,t] = ts_sol[t]
                end
            end
        end
    end
    return nothing
end


function corr_forcing_bias!(pred_frc_corr,obs_t,obs_frc,pred_t,pred_frc)
    # =========================================================================================== 
    # Bias correction of each forcing term
    # The mean seasonal cycle of the ECCO and predicted forcing is estimated,
    # the predicted mean seasonal cycle is removed and replaced by the ECCO mean seasonal cycle
    # That means that the predicted forcing contains the ECCO MSC and anomalies w.r.t. the cycle
    # ===========================================================================================
    println("    Bias correction...")
    # Time steps
    obs_month      = zeros(Int32,length(obs_t))
    pred_month     = zeros(Int32,length(pred_t))
    [obs_month[i]  = Month(obs_t[i]).value for i ∈ 1:length(obs_t)]
    [pred_month[i] = Month(pred_t[i]).value for i ∈ 1:length(pred_t)]
    ens_size = size(pred_frc,4)
    Threads.@threads for mnth ∈ 1:12
        pred_frc_lcl = @views pred_frc[:,:,:,:,pred_month.==mnth];
        pred_frc_corr_lcl = @views pred_frc_corr[:,:,:,:,pred_month.==mnth];
        obs_frc_lcl  = @views obs_frc[:,:,:,obs_month.==mnth];
        @fastmath @inbounds for k ∈ 1:13
            for j ∈ 1:90
                for i ∈ 1:90
                    obs_mean = mean(obs_frc_lcl[i,j,k,:])
                    for ens ∈ 1:ens_size
                        pred_mean = mean(pred_frc_lcl[i,j,k,ens,:])
                        for tstep ∈ 1:size(pred_frc_lcl,5)
                            pred_frc_corr_lcl[i,j,k,ens,tstep] = pred_frc_lcl[i,j,k,ens,tstep] - pred_mean + obs_mean
                        end
                    end
                end
            end
        end
    end
    return nothing
end

function corr_forcing_bias_no_ovl!(pred_frc_corr,obs_t,obs_frc,pred_t,pred_frc,tsteps_obs_inc)
    # =========================================================================================== 
    # Bias correction of each forcing term
    # The mean seasonal cycle of the ECCO and predicted forcing is estimated,
    # the predicted mean seasonal cycle is removed and replaced by the ECCO mean seasonal cycle
    # That means that the predicted forcing contains the ECCO MSC and anomalies w.r.t. the cycle
    # ===========================================================================================
    # Time steps
    obs_month      = zeros(Int32,length(obs_t))
    pred_month     = zeros(Int32,length(pred_t))
    [obs_month[i]  = Month(obs_t[i]).value for i ∈ 1:length(obs_t)]
    [pred_month[i] = Month(pred_t[i]).value for i ∈ 1:length(pred_t)]
    @. obs_month[~tsteps_obs_inc] = 99
    ens_size = size(pred_frc,4)
    Threads.@threads for mnth in 1:12
        pred_frc_lcl = @views pred_frc[:,:,:,:,pred_month.==mnth];
        pred_frc_corr_lcl = @views pred_frc_corr[:,:,:,:,pred_month.==mnth];
        obs_frc_lcl  = @views obs_frc[:,:,:,obs_month.==mnth];
        @fastmath @inbounds for k ∈ 1:13
            for j in 1:90
                for i in 1:90
                    obs_mean = mean(obs_frc_lcl[i,j,k,:])
                    for ens in 1:ens_size
                        pred_mean = mean(pred_frc_lcl[i,j,k,ens,:])
                        for tstep in 1:size(pred_frc_lcl,5)
                            pred_frc_corr_lcl[i,j,k,ens,tstep] = pred_frc_lcl[i,j,k,ens,tstep] - pred_mean + obs_mean
                        end
                    end
                end
            end
        end
    end
    return nothing
end

function read_forcing_file!(obs_frc, idx, filename)
    # ------------------------------------------------------------
    # Read a LLC90 binary forcing file and un-spaghetti tiles 8-13
    # Input:
    #  ̇ filename: string of the filename to read
    # Output:
    # ⋅ store in 'forcing' array at time index 'idx'
    # ------------------------------------------------------------
    n_vars = stat(filename).size
    @assert (n_vars == 4*13*90*90)

    llc_raw = zeros(Float32,(90,90*13));
    read!(filename,llc_raw);
    @. llc_raw = ntoh(llc_raw);    

    @inbounds for n in 1:7
        nstart = 90*(n-1) + 1
        nstop  = nstart+90-1
        @. obs_frc[:,:,n,idx] = llc_raw[:,nstart:nstop]
    end
    @inbounds for n in 8:10
        idx_n = [1:3:270...] .+ (7*90+n-8)
        @. obs_frc[:,:,n,idx] = llc_raw[:,idx_n]
    end
    @inbounds for n in 11:13
        idx_n = [1:3:270...] .+ (10*90+n-11)
        @. obs_frc[:,:,n,idx] = llc_raw[:,idx_n]
    end
    return nothing
end

struct Sensitivity
    time::Array{Date,1}
    sens::Array{Float32,4}
end

function read_sensitivity!(sensitivity,proc_name,settings)
    println("    Reading adjoint sensitivity...")
    dir_sensitivities = homedir()*"/Data/ECCO/RISE/adjoint_sensitivities/"*settings["location"]*"/"
    yr=2004
    Threads.@threads for mnth in 1:12
        filename = dir_sensitivities*lpad(string(mnth),2,"0")*lpad(string(yr-2000),2,"0")*"/adxx_"*proc_name*".0000000129_dim.data";
        sens_lcl = read_sensitivity_file(filename);
        ts_final = sens_ts_final(mnth)
        sens_ts = ts_final .- [size(sens_lcl,4)-1:-1:0...] .* Week(1) .- Day(3)
        sens_smooth = similar(sens_lcl)
        sens_smooth[:,:,:,1] = @views sens_lcl[:,:,:,1];
        for t=2:size(sens_smooth,4)
            @. sens_smooth[:,:,:,t] = 0.5*(@views sens_lcl[:,:,:,t-1] + @views sens_lcl[:,:,:,t])
        end
        sensitivity[mnth] = Sensitivity(sens_ts,sens_smooth)
    end
    return nothing
end

function read_sensitivity_file(filename)
    # ----------------------------------------------------
    # Read a LLC90 binary file and un-spaghetti tiles 8-13
    # Input:
    #  ̇ filename: string of the filename to read
    # Output:
    # ⋅ sensitivity(90,90,13,number of time steps): data
    # ---------------------------------------------------- 
    n_tsteps = convert(Int,stat(filename).size/4/13/90/90)
    llc_raw = zeros(Float32,(90,90*13,n_tsteps));
    read!(filename,llc_raw);
    @. llc_raw = ntoh(llc_raw);    

    sensitivity = zeros(Float32,(90,90,13,n_tsteps));
    @inbounds for n in 1:7
        nstart = 90*(n-1) + 1
        nstop  = nstart+90-1
        sensitivity[:,:,n,:] = @views llc_raw[:,nstart:nstop,:]
    end
    for n in 8:10
        idx_n = [1:3:270...] .+ (7*90+n-8)
        sensitivity[:,:,n,:] = @views llc_raw[:,idx_n,:]
    end

    for n in 11:13
        idx_n = [1:3:270...] .+ (10*90+n-11)
        sensitivity[:,:,n,:] = @views llc_raw[:,idx_n,:]
    end
    return sensitivity
end

function sens_ts_final(mnth)
    if mnth==1
        ts_final = DateTime(2004,2,3)
    elseif mnth==2
        ts_final = DateTime(2004,3,2)
    elseif mnth==3
        ts_final = DateTime(2004,4,6)
    elseif mnth==4
        ts_final = DateTime(2004,5,4)
    elseif mnth==5
        ts_final = DateTime(2004,6,1)
    elseif mnth==6
        ts_final = DateTime(2004,7,6)
    elseif mnth==7
        ts_final = DateTime(2004,8,3)
    elseif mnth==8
        ts_final = DateTime(2004,9,7)
    elseif mnth==9
        ts_final = DateTime(2004,10,5)
    elseif mnth==10
        ts_final = DateTime(2004,11,2)
    elseif mnth==11
        ts_final = DateTime(2004,12,7)
    elseif mnth==12
        ts_final = DateTime(2005,1,4)
    end
    return ts_final
end

function def_sens_timesteps()
    sens_timesteps = zeros(Date,12,2)
    # Jan 2004
    sens_timesteps[1,1] = DateTime(2003,1,29)
    sens_timesteps[1,2] = DateTime(2004,2,4)
    # Feb 2004
    sens_timesteps[2,1] = DateTime(2003,2,26)
    sens_timesteps[2,2] = DateTime(2004,3,3)
    # Mar 2004
    sens_timesteps[3,1] = DateTime(2003,3,26)
    sens_timesteps[3,2] = DateTime(2004,4,7)
    # Apr 2004
    sens_timesteps[4,1] = DateTime(2003,4,30)
    sens_timesteps[4,2] = DateTime(2004,5,5)
    # May 2004
    sens_timesteps[5,1] = DateTime(2003,5,28)
    sens_timesteps[5,2] = DateTime(2004,6,2)
    # June 2004
    sens_timesteps[6,1] = DateTime(2003,6,25)
    sens_timesteps[6,2] = DateTime(2004,7,7)
    # July 2004
    sens_timesteps[7,1] = DateTime(2003,7,30)
    sens_timesteps[7,2] = DateTime(2004,8,4)
    # Aug 2004
    sens_timesteps[8,1] = DateTime(2003,8,27)
    sens_timesteps[8,2] = DateTime(2004,9,1)
    # Sep 2004
    sens_timesteps[9,1] = DateTime(2003,10,1)
    sens_timesteps[9,2] = DateTime(2004,10,6)
    # Oct 2004
    sens_timesteps[10,1] = DateTime(2003,10,29)
    sens_timesteps[10,2] = DateTime(2004,11,3)
    # Nov 2004
    sens_timesteps[11,1] = DateTime(2003,11,26)
    sens_timesteps[11,2] = DateTime(2004,12,1)
    # Dec 2004
    sens_timesteps[12,1] = DateTime(2003,12,31)
    sens_timesteps[12,2] = DateTime(2005,1,5)
    return sens_timesteps
end

function save_reconstruction(conv_obs,conv_seas,conv_obs_seas,conv_obs_pred,ic_mnth,settings)
    println("Saving...")
    # Save reconstruction into netCDF format
    if settings["exclude_pred_mnth_in_bc"]
        fn_save = homedir()*"/Data/ECCO/RISE/results/"*settings["location"]*"/"*settings["location"]*"_pred_"*settings["pred_model"]*"_excl_pred_"*settings["months"][ic_mnth]*"IC.nc"        
    else
        fn_save = homedir()*"/Data/ECCO/RISE/results/"*settings["location"]*"/"*settings["location"]*"_pred_"*settings["pred_model"]*"_"*settings["months"][ic_mnth]*"IC.nc"        
    end
    fh = Dataset(fn_save,"c")
    defDim(fh,"year",length(settings["years"]))
    defDim(fh,"time",size(settings["time"],2))
    defDim(fh,"month",12)
    defDim(fh,"forcing",4)
    defDim(fh,"ens_member",length(settings["ens_range"]))
    defVar(fh,"year",settings["years"] ,("year",),deflatelevel=4)
    defVar(fh,"time", convert.(DateTime,settings["time"][ic_mnth,:]),("time",),deflatelevel=4)
    defVar(fh,"month",[1:12...] ,("month",),deflatelevel=4)
    defVar(fh,"forcing",settings["processes"] ,("forcing",),deflatelevel=4)
    defVar(fh,"ens_member",settings["ens_range"],("ens_member",),deflatelevel=4,attrib = Dict("long_name" => "Ensemble member"))
    defVar(fh,"time_init", convert.(DateTime,unique(settings["time_start"][ic_mnth,:])),("year",),deflatelevel=4,attrib = Dict("long_name" => "Time step at which CCSM run is initialized"))
    defVar(fh,"conv_obs",conv_obs,("forcing","time"),fillvalue=-9999,deflatelevel=4,attrib = Dict("units" => "m","long_name" => "Reconstruction from convolving ECCO forcing for all timesteps"))
    defVar(fh,"conv_seas",conv_seas,("forcing","time"),fillvalue=-9999,deflatelevel=4,attrib = Dict("units" => "m","long_name" => "Forecast from convolving ECCO mean seasonal cycle forcing for all time steps"))
    defVar(fh,"conv_obs_seas",conv_obs_seas,("forcing","time"),fillvalue=-9999,deflatelevel=4,attrib = Dict("units" => "m","long_name" => "Forecast from convolving ECCO forcing for timesteps before initialization time and ECCO mean seasonal cycle forcing for timesteps after initialization time"))
    defVar(fh,"conv_obs_pred",conv_obs_pred,("forcing","ens_member","time"),fillvalue=-9999,deflatelevel=4,attrib = Dict("units" => "m","long_name" => "Forecast from convolving ECCO forcing for timesteps before initialization time and CCSM forcing for timesteps after initialization time"))
    close(fh)
    return nothing
end

main()