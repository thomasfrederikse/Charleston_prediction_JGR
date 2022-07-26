# Read observed sea level in Charleston and make projections using autocorrelation/persistence
using StatsBase
using Statistics
using DelimitedFiles
using NCDatasets
using NetCDF
using Interpolations
using Dates
function main()
    settings = def_settings()
    obs_time,obs_msl = read_obs(settings)# Read tide-gauge data and remove IB
    tvec = zeros(Int,length(obs_time),2)
    @. tvec[:,1] = floor(obs_time)
    @. tvec[:,2] = @. round(12*(obs_time - floor(obs_time) + 1/24 ))

    # Remove trend and seasonal cycle
    amat = ones(length(obs_time),8)
    amat[:,2] = obs_time .- mean(obs_time)
    @. amat[:,3] = sin(2*π*obs_time)
    @. amat[:,4] = cos(2*π*obs_time)
    @. amat[:,5] = sin(4*π*obs_time)
    @. amat[:,6] = cos(4*π*obs_time)
    @. amat[:,7] = sin(6*π*obs_time)
    @. amat[:,8] = cos(6*π*obs_time)
    sol = amat\obs_msl
    msl_corr = obs_msl - amat*sol
    
    pred_array = zeros(12,12*length(settings["years"]))
    for ic_mnth ∈ 1:12
        ts = 12 + ic_mnth
        tsteps = [ts:12:length(tval)...]
        ϕ₁ = msl_corr[tsteps.-1]\msl_corr[tsteps]
        println(ϕ₁)
        for (yr_idx,yr) ∈ enumerate(settings["years"])
            proj_yr = zeros(12)
            start_idx = findfirst((tvec[:,1] .== yr ) .& (tvec[:,2] .== ic_mnth))
            proj_yr[1] = ϕ₁ * msl_corr[start_idx-1]
            for lag ∈ 2:12
                proj_yr[lag] = ϕ₁ * proj_yr[lag-1]
            end
            pred_array[ic_mnth,12*(yr_idx-1)+1:12*yr_idx] = proj_yr
        end
    end

    # Save predictions
    tvec_save = zeros(Int,12,12*length(settings["years"]),2)
    for ic_mnth ∈ 1:12
        tdata = [Date(settings["years"][1],ic_mnth,15):Month(1):Date(2022,1,15)...][1:(12*length(settings["years"]))]
        tvec_save[ic_mnth,:,1] = [Year(tdata[i]).value for i ∈ 1:12*length(settings["years"])]
        tvec_save[ic_mnth,:,2] = [Month(tdata[i]).value for i ∈ 1:12*length(settings["years"])]
    end

    fh = Dataset(settings["fn_persistence"],"c")
    defDim(fh,"IC_mnth",12)
    defDim(fh,"t",size(tvec_save,2))
    defDim(fh,"dbl",2)

    defVar(fh,"IC_mnth",[1:12...],("IC_mnth",),deflatelevel=4)
    defVar(fh,"tvec",tvec_save,("tstart","t","dbl"),deflatelevel=4)

    defVar(fh,"eta_persistence",pred_array,("IC_mnth","t"),deflatelevel=4)
    close(fh)
end

function def_settings()
    settings=Dict()
    settings["location"] = "SanDiego"
    settings["years"]  = [1995:2016...]
    settings["mode"] = "TG"

    # Directories
    settings["dir_data"] = homedir()*"/Data/ECCO/RISE/"
    settings["dir_results"] = settings["dir_data"]*"results/"*settings["location"]*"/"
    settings["fn_persistence"] = settings["dir_results"]*"/"*settings["location"]*"_persistence.nc"
    settings["fn_observations"] = settings["dir_results"]*settings["location"]*"_observations.nc"
    return settings
end

function read_obs(settings)
    obs_time = ncread(settings["fn_observations"],"time")
    obs_msl = ncread(settings["fn_observations"],"TG_lmsl")
    return obs_time,obs_msl
end
