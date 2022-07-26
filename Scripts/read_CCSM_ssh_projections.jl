# --------------------------------------------------------------------------
# Read the CCSM4 SSH projections for each ensemble member and starting month
# - Remove drift
# - Save as netCDF
# --------------------------------------------------------------------------
using Plots
using NetCDF
using Statistics
using DelimitedFiles
using NCDatasets
using Dates
function main()
    # Settings
    settings = Dict()
    settings["location"] = "Charleston"
    # Folders
    settings["dir_ccsm"] = homedir()*"/Data/ECCO/RISE/prediction_models/ccsm/raw/"
    settings["fn_ccsm_ssh"] = homedir()*"/Data/ECCO/RISE/results/"*settings["location"]*"/"*settings["location"]*"_ccsm_ssh.nc"
    settings["fn_observations"] = homedir()*"/Data/ECCO/RISE/results/"*settings["location"]*"/"*settings["location"]*"_observations.nc"
    # Time steps and ensembles
    settings["months"]   = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
    settings["ens_range"]  = [1:10...] # Array that sets which ensemble members to include
    settings["years"] = [1994:2016...]
    # Coordinates from location file
    settings["ϕ"] = ncread(settings["fn_observations"],"lon")[1]
    settings["θ"] = ncread(settings["fn_observations"],"lat")[1]
    
    tvec,η_full = read_CCSM_ssh_raw(settings)
    η_dedrift = remove_seasonal_drift(tvec,η_full,settings)
    save_data(tvec,η_full,η_dedrift,settings)
    return nothing
end

function read_CCSM_ssh_raw(settings)
    η_full = zeros(length(settings["months"]),length(settings["ens_range"]),12*length(settings["years"]))
    # Vector with time step info
    tvec = zeros(Int,12,12*length(settings["years"]),2)
    for ic_mnth ∈ 1:length(settings["months"])
        tdata = [Date(settings["years"][1],ic_mnth,15):Month(1):Date(2022,1,15)...][1:(12*length(settings["years"]))]
        tvec[ic_mnth,:,1] = [Year(tdata[i]).value for i ∈ 1:12*length(settings["years"])]
        tvec[ic_mnth,:,2] = [Month(tdata[i]).value for i ∈ 1:12*length(settings["years"])]
    end

    # Select grid cell to read
    fn = settings["dir_ccsm"] * settings["months"][1] *"IC/ccsm4_0_cfsrr_Fcst.E"*string(1)*".pop.h."*string(tvec[1,1,1])*"-"*lpad(tvec[1,1,2],2,"0")*".nc"
    ϕ=ncread(fn,"lon")    
    θ=ncread(fn,"lat")
    slm = ncread(fn,"SSH")[:,:,1] .< 1e30
    dst = sqrt.(reshape((settings["ϕ"].-ϕ).^2,(:,1)) .+ reshape((settings["θ"].-θ).^2,(1,:)))
    @. dst[!slm] = 1e10
    locs = argmin(dst)
    ϕᵢ = locs[1]
    θᵢ = locs[2]

    for ic_mnth ∈ 1:length(settings["months"])
        for ens ∈ settings["ens_range"]
            for tstep ∈ 1:12*length(settings["years"])
                fn = settings["dir_ccsm"] * settings["months"][ic_mnth] *"IC/ccsm4_0_cfsrr_Fcst.E"*string(ens)*".pop.h."*string(tvec[ic_mnth,tstep,1])*"-"*lpad(tvec[ic_mnth,tstep,2],2,"0")*".nc"
                η_full[ic_mnth,ens,tstep] = ncread(fn,"SSH",start=[ϕᵢ,θᵢ,1],count=[1,1,1])[1] ./100 # From cm to m
            end
        end
    end

    # Remove mean
    for ic_mnth ∈ 1:length(settings["months"])
        for ens ∈ settings["ens_range"]
            η_full[ic_mnth,ens,:] .-= mean(η_full[ic_mnth,ens,:])
        end
    end
    return tvec,η_full
end

function remove_seasonal_drift(tvec,η_full,settings)
    η_dedrift = similar(η_full)
    for ic_mnth ∈ 1:length(settings["months"])
        for ts_mnth ∈ 1:12
            mnth_idx = tvec[ic_mnth,:,2] .== ts_mnth
            for ens ∈ settings["ens_range"]
                η_dedrift[ic_mnth,ens,mnth_idx] = η_full[ic_mnth,ens,mnth_idx] .- mean(η_full[ic_mnth,ens,mnth_idx])
            end
        end
    end
    return η_dedrift
end

function location_list(loc)
    if loc=="Charleston"
        θ = 32.782
        ϕ = 360-79.925
    elseif loc=="SanDiego"
        θ=32.866667
        ϕ=360-117.256667
    else
        error("Location "*loc*" not recognized")
    end
    return ϕ,θ
end

function save_data(tvec,η_full,η_dedrift,settings)
    fh = Dataset(settings["fn_ccsm_ssh"],"c")
    defDim(fh,"IC_mnth",length(settings["months"]))
    defDim(fh,"ens",length(settings["ens_range"]))
    defDim(fh,"t",size(tvec,2))
    defDim(fh,"dbl",2)

    defVar(fh,"ens",settings["ens_range"],("ens",),deflatelevel=4)
    defVar(fh,"IC_mnth",[1:12...],("IC_mnth",),deflatelevel=4)
    defVar(fh,"tvec",tvec,("tstart","t","dbl"),deflatelevel=4)

    defVar(fh,"eta_full",η_full,("IC_mnth","ens","t"),deflatelevel=4)
    defVar(fh,"eta_dedrift",η_dedrift,("IC_mnth","ens","t"),deflatelevel=4)
    close(fh)
    return nothing
end

