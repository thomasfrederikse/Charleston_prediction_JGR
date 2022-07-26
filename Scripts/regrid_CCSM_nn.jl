# --------------------------------------------
# Regrid CCSM data onto the ECCO LLC90 grid
# Use nearest-neighboor interpolation
# Save in same format as CCSM regridded netCDF
# --------------------------------------------
using Distributed
@everywhere using NCDatasets
using NetCDF
using NearestNeighbors
using Dates

include("regrid_functions.jl")

function main()
    settings = DefSettings() 
    RegridCCSM(settings)
    return nothing
end

function DefSettings()
    settings = Dict()
    settings["model_name"] = "ccsm"
    settings["process_ecco"]     = ["empmr","qnet","tauu","tauv"]
    settings["conversion_factor"] = Dict("empmr"=>-1.0f0, "qnet"=>-1.0f0, "tauu"=>-0.1f0, "tauv"=>-0.1f0, )

    settings["init_month_list"] = [1:12...]
    settings["ens_list"]   = [1:10...]
    settings["init_months"] = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"][settings["init_month_list"]]

    settings["dir_raw"]    = homedir()*"/Data/ECCO/RISE/prediction_models/ccsm/raw/"
    settings["dir_regrid"] = homedir()*"/Data/ECCO/RISE/prediction_models/ccsm/regrid_nn/"

    settings["dir_ECCO"] =  homedir()*"/Data/ECCO/v4r4/";
    settings["fn_ECCO_grid"] = settings["dir_ECCO"]*"ECCOv4r4_grid.nc"

    settings["ylist"] = [1994:2020...]
    return settings
end

function RegridCCSM(settings)
    # Prepare NN trees
    acc_ccsm, lon_ccsm, lat_ccsm = read_ccsm_coords(settings)
    kdtree_ccsm = KDTree([lon_ccsm lat_ccsm]',leafsize=1);

    # Read ECCO grid
    lon_ECCO,lat_ECCO,slm_ECCO = read_ECCO_grid(settings);
    
    # Process start times and ensembles
    for init_month ∈ settings["init_month_list"]
        for ens ∈ settings["ens_list"]
            println(" Init month "*string(init_month)*" - ensemble "*string(ens)*"...")
            tvec,empmr_ccsm,qnet_ccsm,tauu_ccsm,tauv_ccsm = read_ccsm_data(init_month,ens,acc_ccsm,settings)
            
            # Pre-allocate re-gridded data
            empmr_ECCO = zeros(Float32,size(lon_ECCO,1),size(lon_ECCO,2),size(lon_ECCO,3),size(tvec,1));
            qnet_ECCO = zeros(Float32,size(lon_ECCO,1),size(lon_ECCO,2),size(lon_ECCO,3),size(tvec,1));
            tauu_ECCO = zeros(Float32,size(lon_ECCO,1),size(lon_ECCO,2),size(lon_ECCO,3),size(tvec,1));
            tauv_ECCO = zeros(Float32,size(lon_ECCO,1),size(lon_ECCO,2),size(lon_ECCO,3),size(tvec,1));
            
            # Interpolate
            Threads.@threads for k ∈ 1:size(lon_ECCO,3)
                for j ∈ 1:size(lon_ECCO,2)
                    for i ∈ 1:size(lon_ECCO,1)
                        idx_reg = nn(kdtree_ccsm,[lon_ECCO[i,j,k],lat_ECCO[i,j,k]])[1]
                        empmr_ECCO[i,j,k,:] = @views empmr_ccsm[idx_reg,:]
                        qnet_ECCO[i,j,k,:]  = @views qnet_ccsm[idx_reg,:]
                        tauu_ECCO[i,j,k,:]  = @views tauu_ccsm[idx_reg,:]
                        tauv_ECCO[i,j,k,:]  = @views tauv_ccsm[idx_reg,:]
                    end
                end
            end

            # Mask out land
            @. empmr_ECCO[!slm_ECCO,:] = 0;
            @. qnet_ECCO[!slm_ECCO,:] = 0;
            @. tauu_ECCO[!slm_ECCO,:] = 0;
            @. tauv_ECCO[!slm_ECCO,:] = 0;

            # Save
            @spawnat :any save_file(lon_ECCO,lat_ECCO,empmr_ECCO,qnet_ECCO,tauu_ECCO,tauv_ECCO,init_month,ens,tvec,settings)
        end
    end
    return 
end

function read_ccsm_coords(settings)
    fn = settings["dir_raw"]*"/JanIC/ccsm4_0_cfsrr_Fcst.E4.pop.h.2020-12.nc"
    acc_ccsm = ncread(fn,"TAUX")[:,:,1][:] .< 1e30
    lon_ccsm = repeat(ncread(fn,"lon"),1,180)[acc_ccsm];
    lat_ccsm = repeat(reshape(ncread(fn,"lat"),(1,:)),360,1)[acc_ccsm];
    return acc_ccsm, lon_ccsm, lat_ccsm
end

function read_ccsm_data(init_month,ens,acc_ccsm,settings)
    # Define time steps
    tvec = zeros(Int32,12*length(settings["ylist"]),2)
    tvec[1,1] = settings["ylist"][1]
    tvec[:,2] = [1:size(tvec,1)...] .+ (init_month-1)
    tdiv = divrem.(tvec[:,2].-1,12)
    for i ∈ 1:size(tvec,1)
        tvec[i,1] = settings["ylist"][1] + tdiv[i][1]
        tvec[i,2] = tdiv[i][2]+1
    end

    # Allocate arrays
    empmr_ccsm = zeros(Float32,sum(acc_ccsm),size(tvec,1))
    qnet_ccsm = zeros(Float32,sum(acc_ccsm),size(tvec,1))
    tauu_ccsm = zeros(Float32,sum(acc_ccsm),size(tvec,1))
    tauv_ccsm = zeros(Float32,sum(acc_ccsm),size(tvec,1))

    for tstep ∈ 1:size(tvec,1)
        fn = settings["dir_raw"]*settings["init_months"][init_month]*"IC/ccsm4_0_cfsrr_Fcst.E"*string(ens)*".pop.h."*string(tvec[tstep,1])*"-"*lpad(tvec[tstep,2],2,"0")*".nc"
        empmr_ccsm[:,tstep] = settings["conversion_factor"]["empmr"] .* (ncread(fn,"EVAP_F")[acc_ccsm] .+ ncread(fn,"PREC_F")[acc_ccsm]) 
        qnet_ccsm[:,tstep] = settings["conversion_factor"]["qnet"] .* ncread(fn,"SHF")[acc_ccsm] 
        tauu_ccsm[:,tstep] = settings["conversion_factor"]["tauu"] .* ncread(fn,"TAUX")[acc_ccsm] 
        tauv_ccsm[:,tstep] = settings["conversion_factor"]["tauv"] .* ncread(fn,"TAUY")[acc_ccsm] 
    end
    return(tvec,empmr_ccsm,qnet_ccsm,tauu_ccsm,tauv_ccsm)
end

main()