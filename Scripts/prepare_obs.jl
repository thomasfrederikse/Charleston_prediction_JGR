# -----------------------------------
# Prepare observations for comparison
# in seasonal predictions
# -----------------------------------
# 1. Tide gauge data
# 2. Inverse-barometer correction
# 3. Altimetry
# 4. Estimate altimetry GMSL for 
#    removal from TG and altimetry
#    observations
# -----------------------------------
using NetCDF
using NCDatasets
using Statistics
using Interpolations
using Dates
using DelimitedFiles

function main()
    settings = def_settings()
    AVISO = read_AVISO(settings)
    TG = read_TG(settings)

    # Remove GMSL over altimetry period
    TG["lmsl"] = TG["lmsl"] - AVISO["gmsl"] - TG["IB"]
    AVISO["lmsl"] = AVISO["lmsl"] - AVISO["gmsl"]
    TG["lmsl_extend"] = TG["lmsl_extend"] - TG["IB_extend"]
    save_data(AVISO,TG,settings)
    return nothing
end

function def_settings()
    settings=Dict()
    settings["location"] = "SanDiego"
    settings["ϕ"], settings["θ"], settings["PSMSL_id"] = location_list(settings["location"])
    settings["time"] = [1993.0f0+1f0/24f0:1f0/12f0:2021f0-1f0/24f0...]
    settings["time_extend"] = [1991.0f0+1f0/24f0:1f0/12f0:2021f0-1f0/24f0...]
    settings["dir_data"] = homedir()*"/Data/ECCO/RISE/"
    settings["fn_save"] = settings["dir_data"]*"results/"*settings["location"]*"/"*settings["location"]*"_observations.nc"
    return settings
end

function location_list(loc)
    if loc=="Charleston"
        θ = 32.782
        ϕ = 360-79.925
        PSMSL_id = 234
    elseif loc=="SanDiego"
        θ=32.866667
        ϕ=360-117.256667
        PSMSL_id = 158
    else
        error("Location "*loc*" not recognized")
    end
    return ϕ,θ,PSMSL_id
end

function read_AVISO(settings)
    AVISO = Dict()
    fn = homedir()*"/Data/Altimetry/CDS/CDS_monthly_1993_2020.nc"
    ϕ = ncread(fn,"lon")
    θ = ncread(fn,"lat")
    time = ncread(fn,"time")

    ϕ_idx = argmin(@.(abs(settings["ϕ"] - ϕ)))
    θ_idx = argmin(@.(abs(settings["θ"] - θ)))

    AVISO["lmsl"]  = LinearInterpolation(time,ncread(fn,"ssh",start=[ϕ_idx,θ_idx,1],count=[1,1,-1])[1,1,:] .* 0.05,extrapolation_bc = NaN)(settings["time"]) ./1000 # Units m
    AVISO["gmsl"] = LinearInterpolation(time,ncread(fn,"gmsl"),extrapolation_bc = NaN)(settings["time"]) ./1000 # Units m

    AVISO["lmsl"] .-= mean(filter(isfinite,AVISO["lmsl"]))
    AVISO["gmsl"] .-= mean(filter(isfinite,AVISO["gmsl"]))
    return AVISO
end

function read_TG(settings)
    η_IB,η_IB_extend = compute_IB(settings)
    fn_psmsl = homedir()*"/Data/TideGauges/rlr_monthly/data/"*string(settings["PSMSL_id"])*".rlrdata"
    TG_raw = readdlm(fn_psmsl,';')[:,1:2]
    acc_idx = TG_raw[:,2] .> 0
    TG_int = LinearInterpolation(TG_raw[acc_idx,1],TG_raw[acc_idx,2])(settings["time"])
    TG_int = (TG_int .- mean(TG_int)) ./ 1000
    TG_int_extend = LinearInterpolation(TG_raw[acc_idx,1],TG_raw[acc_idx,2])(settings["time_extend"])
    TG_int_extend = (TG_int_extend .- mean(TG_int_extend)) ./ 1000

    TG = Dict()
    TG["lmsl"] = TG_int
    TG["lmsl_extend"] = TG_int_extend
    TG["IB"] = η_IB
    TG["IB_extend"] = η_IB_extend
    return TG
end

function compute_IB(settings)
    # Read data
    fn_ERA = homedir()*"/Data/Reanalyses/ERA5/ERA5.nc"
    fh_ERA = Dataset(fn_ERA,"r")
    ϕ = fh_ERA["lon"][:]
    θ = fh_ERA["lat"][:]
    slm = convert.(Bool,fh_ERA["slm"][:])
    close(fh_ERA)

    dst = @. sqrt((settings["ϕ"] - $reshape(ϕ,(:,1)))^2+(settings["θ"] - $reshape(θ,(1,:)))^2)
    @. dst[~slm] = 1.0e5
    loc_idx = argmin(dst)

    fh_ERA = Dataset(fn_ERA,"r")
    time = fh_ERA["time"][:]
    acc_time = time .< DateTime(2021,1,1)
    time = time[acc_time]
    mslp = nomissing(fh_ERA["msl"][loc_idx[1],loc_idx[2],:],NaN)[acc_time]
    mslp_ocean = nomissing(fh_ERA["mslp_ocean"][:],NaN)[acc_time]
    tval = [Year(time[t]).value+Month(time[t]).value/12-1/24 for t ∈ 1:length(time)]

    ρ = 1002 # Sea water density
    g = 9.79 # Earth gravity
    η_IB_raw = (mslp .- mslp_ocean) ./ (-(ρ * g)); # IB in meters
    η_IB_raw .-= mean(η_IB_raw)
    η_IB = LinearInterpolation(tval,η_IB_raw)(settings["time"])
    η_IB_extend = LinearInterpolation(tval,η_IB_raw)(settings["time_extend"])
    return η_IB,η_IB_extend
end


function grid_area(lon,lat)
    gridsize = abs(lat[2]-lat[1])
    area = @. deg2rad(gridsize) * (sind(lat+gridsize/2)-sind(lat-gridsize/2)) * 6371000^2
    return repeat(area',size(lon,1))
end

function save_data(AVISO,TG,settings)
    fh = Dataset(settings["fn_save"],"c")
    defDim(fh,"time",length(settings["time"]))
    defDim(fh,"time_long",length(settings["time_extend"]))
    defDim(fh,"coords",1)

    defVar(fh,"time",settings["time"],("time",),deflatelevel=4)
    defVar(fh,"time_long",settings["time_extend"],("time_long",),deflatelevel=4)
    defVar(fh,"lon",[settings["ϕ"]],("coords",),deflatelevel=4)
    defVar(fh,"lat",[settings["θ"]],("coords",),deflatelevel=4)
    defVar(fh,"location",[settings["location"]],("coords",),deflatelevel=4)

    defVar(fh,"Alt_lmsl",AVISO["lmsl"],("time",),deflatelevel=4)
    defVar(fh,"TG_lmsl",TG["lmsl"],("time",),deflatelevel=4)
    defVar(fh,"gmsl",AVISO["gmsl"],("time",),deflatelevel=4)
    defVar(fh,"IB",TG["IB"],("time",),deflatelevel=4)
    defVar(fh,"IB_long",TG["IB_extend"],("time_long",),deflatelevel=4)
    defVar(fh,"TG_lmsl_long",TG["lmsl_extend"],("time_long",),deflatelevel=4)
    close(fh)
    return nothing
end

function junk()
    scatter([1994:2015...],ccsm_mn[1:12:end],marker=:circ,label=false,size=(1000,500),margin=5Plots.mm)
    plot!(settings["time"],lmsl_noseas.-IB_noseas,linewidth=2,label="Tide gauge obs")
    plot!(t,mean(η_dedrift[1,:,:],dims=1)[1,:],linewidth=2,label="CCSM4 dedrifted")
    ylabel!("Height (m)")
    savefig("CCSM4_1994_2015_JanIC.png")
end