# -----------------------------------------------------
# Store generic regridding functions that can be shared
# among different data sets
# -----------------------------------------------------
@everywhere function save_file(lon_ECCO,lat_ECCO,empmr_ECCO,qnet_ECCO,tauu_ECCO,tauv_ECCO,init_month,ens,tvec,settings)
    fn_save = settings["dir_regrid"]*settings["model_name"]*"_init_"*lpad(init_month,2,"0")*"_ens_"*lpad(ens,2,"0")*".nc"
    fh = Dataset(fn_save,"c")
    t = convert.(Float32,tvec[:,1] .+ (tvec[:,2].-0.5) ./12)
    defDim(fh,"time",length(t))
    defDim(fh,"yrmn",2)
    defDim(fh,"x",size(lon_ECCO,1))
    defDim(fh,"y",size(lon_ECCO,2))
    defDim(fh,"tile",size(lon_ECCO,3))

    defVar(fh,"time",t,("time",),fillvalue=-9999,deflatelevel=2)
    defVar(fh,"tvec",tvec,("time","yrmn",),fillvalue=-9999,deflatelevel=2)
    defVar(fh,"lon",lon_ECCO,("x","y","tile"),fillvalue=-9999,deflatelevel=2)
    defVar(fh,"lat",lat_ECCO,("x","y","tile"),fillvalue=-9999,deflatelevel=2)
    defVar(fh,"empmr",empmr_ECCO,("x","y","tile","time",),deflatelevel=2) 
    defVar(fh,"qnet",qnet_ECCO,("x","y","tile","time",),deflatelevel=2) 
    defVar(fh,"tauu",tauu_ECCO,("x","y","tile","time",),deflatelevel=2) 
    defVar(fh,"tauv",tauv_ECCO,("x","y","tile","time",),deflatelevel=2) 
    close(fh)
    return nothing
end

@everywhere function read_ECCO_grid(settings)
    lon_ECCO = mod.(ncread(settings["fn_ECCO_grid"],"XC"),360)
    lat_ECCO = ncread(settings["fn_ECCO_grid"],"YC")
    slm_ECCO = convert.(Bool,ncread(settings["fn_ECCO_grid"],"hFacC",start=[1,1,1,1],count=[-1,-1,-1,1])[:,:,:,1];)
    return lon_ECCO,lat_ECCO,slm_ECCO
end

@everywhere function prepare_interpolator(lon_pred,lat_pred,lon_ECCO,lat_ECCO)
    scinterp = pyimport("scipy.interpolate")    
    latlon_ECCO=hcat(lon_ECCO[:],lat_ECCO[:])
    itp_ECCO = scinterp.LinearNDInterpolator([lon_pred lat_pred],ones(length(lon_pred)),fill_value=0.0);
    return latlon_ECCO,itp_ECCO
end

@everywhere function interp_on_ecco_grid!(grid_in,itp_function,latlon_ECCO,grid_out)
    for tstep ∈ 1:size(grid_in,2)
        itp_function.values = reshape((@. convert(Float64,@views grid_in[:,tstep])),(:,1));
        grid_out[:,:,:,tstep] = reshape(itp_function(latlon_ECCO),(90,90,13));
    end
    return nothing
end