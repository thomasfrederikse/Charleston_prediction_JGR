# Read bathymetry and SSHDYN
# interpolate onto regular lat lon grid
# Save bathy and mdt for GMT
using Plots
using NCDatasets
using StaticArrays
using LinearAlgebra
using NetCDF
using Dates
using Statistics
using PyCall
using DelimitedFiles

function main()
    tval,ts_ecco_noseas = read_sshdyn_charleston()
    tval,lmsl_psmsl_noseas,lmsl_aviso_noseas = read_obs_data()



    plot(tval,lmsl_psmsl_noseas)
    plot!(tval,lmsl_aviso_noseas)
    plot!(tval,ts_ecco_noseas)

    writedlm(homedir()*"/Projects/2020_ECCO_adjoint/GMT/bathy_mdt_obs/ts_ecco.txt", [tval ts_ecco_noseas], ';');  
    writedlm(homedir()*"/Projects/2020_ECCO_adjoint/GMT/bathy_mdt_obs/ts_aviso.txt", [tval lmsl_aviso_noseas], ';');  
    writedlm(homedir()*"/Projects/2020_ECCO_adjoint/GMT/bathy_mdt_obs/ts_psmsl.txt", [tval lmsl_psmsl_noseas], ';');  
    return nothing
end

function proc_bathy()
    # Read bathymetry
    fn = homedir()*"/Data/ECCO/v4r4/ECCOv4r4_grid.nc"
    depth_llc = ncread(fn,"Depth")
    ϕ,θ,slm,depth = interpolate_2d(depth_llc)
    save_2d(homedir()*"/Projects/2020_ECCO_adjoint/GMT/bathy_mdt_obs/bathy.nc",ϕ,θ,depth[:,:,1])
    return nothing
end

function proc_mdt()
    mdt_llc = zeros(Float32,90,90,13)
    years = [1993:2017...]
    months = [1:12...]
    for yr in years
        for mnth in months
            println(string(yr)*" "*string(mnth))
            fn = homedir()*"/Data/ECCO/v4r4/SSHDYN/"*string(yr)*"/SSHDYN_"*string(yr)*"_"*lpad(string(mnth),2,"0")*".nc"
            mdt_llc += ncread(fn,"SSHDYN");
        end
    end
    mdt_llc ./= (length(years)*length(months))
    ϕ,θ,slm,mdt = interpolate_2d(mdt_llc)
    save_2d(homedir()*"/Projects/2020_ECCO_adjoint/GMT/bathy_mdt_obs/mdt.nc",ϕ,θ,mdt[:,:,1])
    return nothing
end

function save_2d(fn_save,ϕ,θ,fld)
    fh = Dataset(fn_save,"c")
    defDim(fh,"lon",size(ϕ,1))
    defDim(fh,"lat",size(θ,1))
    defVar(fh,"lon",ϕ,("lon",),fillvalue=-9999,deflatelevel=2)
    defVar(fh,"lat",θ,("lat",),fillvalue=-9999,deflatelevel=2)
    defVar(fh,"fld",fld,("lon","lat",),fillvalue=-9999,deflatelevel=2)
    close(fh)
    return nothing
end

function interpolate_2d(field_ecco)
    scinterp = pyimport("scipy.interpolate")    
    # Read ECCO
    fn = homedir()*"/Data/ECCO/v4r4/ECCOv4r4_grid.nc"
    ϕ_ECCO = ncread(fn,"XC");
    θ_ECCO = ncread(fn,"YC");
    slm_ECCO = ncread(fn,"Depth").>0
    ϕ_ECCO[ϕ_ECCO.<0] .+= 360
    ϕ_ECCO = ϕ_ECCO[slm_ECCO]
    θ_ECCO = θ_ECCO[slm_ECCO]
    ϕθ_ECCO=hcat(ϕ_ECCO[:],θ_ECCO[:])
    itp_ecco = scinterp.LinearNDInterpolator(ϕθ_ECCO,ones(size(ϕθ_ECCO,1)),fill_value=0);
    
    ϕ = [0.25:0.5:359.75...]
    θ = [-89.75:0.5:89.75...]
    ϕϕ = zeros(Float32,size(ϕ,1),size(θ,1))
    θθ = zeros(Float32,size(ϕ,1),size(θ,1))
    for i=1:size(ϕ,1),j=1:size(θ,1)
        ϕϕ[i,j] = ϕ[i]
        θθ[i,j] = θ[j]
    end
    ϕθ=hcat(ϕϕ[:],θθ[:])

    field_interp = zeros(size(ϕ,1),size(θ,1));
    for tstep in 1:size(field_ecco,4)
        field_lcl = field_ecco[:,:,:,tstep][slm_ECCO];
        itp_ecco.values = convert.(Float64,reshape(field_lcl,(:,1)));
        field_interp[:,:,tstep] = reshape(itp_ecco(ϕθ),(720,360));
    end

    # Read slm
    slm = 1.0 .- ncread(homedir()*"/Data/GRACE/JPL_mascon/mask.nc","land")
    for tstep in 1:size(field_ecco,4)
        field_interp[:,:,tstep] .*= slm
    end
    field_interp = convert.(Float32,field_interp)
    return(ϕ,θ,slm,field_interp)
end

function read_obs_data()
    fn_aviso = homedir() * "/Projects/2020_ECCO_adjoint/Charleston_TG/AVISO_Charleston_SC.nc"
    fn_psmsl = homedir() * "/Projects/2020_ECCO_adjoint/Charleston_TG/TG_charleston.txt"
    fn_ib    = homedir() * "/Projects/2020_ECCO_adjoint/Charleston_TG/IB_charleston.txt"

    t = Date(1993,1,15) .+ Month.([0:299...])
    tval = [1993+1/24:1/12:2017+23/24...]
    gmsl_aviso = ncread(fn_aviso,"aviso_gmsl")
    lmsl_aviso = ncread(fn_aviso,"aviso_sla_Charleston") .- gmsl_aviso
    IB_charleston = readdlm(fn_ib,';')[169:468,:]    # Load IB
    TG_charleston = readdlm(fn_psmsl,';')[856:end-24,1:2]
    lmsl_psmsl = TG_charleston[:,2] .-IB_charleston[:,3]
    lmsl_psmsl = lmsl_psmsl ./1000 .- gmsl_aviso
    # Separate mean seasonal cycle
    lmsl_psmsl_noseas = separate_mss(t,tval,lmsl_psmsl)
    lmsl_aviso_noseas = separate_mss(t,tval,lmsl_aviso)

    lmsl_psmsl_noseas .-= mean(lmsl_psmsl_noseas)
    lmsl_aviso_noseas .-= mean(lmsl_aviso_noseas)
    return(tval,lmsl_psmsl_noseas,lmsl_aviso_noseas)
end

function separate_mss(t,tval,msl_ts)
    # 1. detrend time series
    amat = ones(size(tval,1),6);        
    amat[:,3] = sin.(2*pi.*tval)
    amat[:,4] = sin.(4*pi.*tval)
    amat[:,5] = cos.(2*pi.*tval)
    amat[:,6] = cos.(4*pi.*tval)    
    amat[:,2] = tval.- mean(tval);
    amat_tr = transpose(amat)
    amat_sq = inv(amat_tr*amat)
    sol = (amat_sq * (amat_tr*msl_ts))
    sol[3:end] .= 0
    msl_ts_detrend = msl_ts - amat * sol

    msl_mnth = zeros(Int,length(t))
    [msl_mnth[i] = Month(t[i]).value for i in 1:length(t)]
    mss = zeros(12)
    mss_ts = zeros(length(t))
    for mnth in 1:12
        month_idx = @. (msl_mnth == mnth)
        mss[mnth] = mean(msl_ts_detrend[month_idx])
        mss_ts[month_idx] .= mss[mnth]
    end
    msl_ts_noseas = msl_ts.-mss_ts
    return msl_ts_noseas
end

function read_sshdyn_charleston()
    years = [1993:2017...]
    months = [1:12...]
    t = Date(1993,1,15) .+ Month.([0:299...])
    tval = [1993+1/24:1/12:2017+23/24...]
    fn = homedir()*"/Data/ECCO/v4r4/ECCOv4r4_grid.nc"

    slm = ((ncread(fn,"hFacS",start=[1,1,1,1],count=[-1,-1,-1,1]) .== 1)  .&  (ncread(fn,"hFacC",start=[1,1,1,1],count=[-1,-1,-1,1]) .== 1) .& (ncread(fn,"hFacW",start=[1,1,1,1],count=[-1,-1,-1,1]) .== 1))[:,:,:,1];
    area = slm .* ncread(fn,"rA")
    θ =ncread(fn,"YC")
    ϕ =ncread(fn,"XC")
    # @. θ[slm] = 1e9;
    # @. ϕ[slm] = 1e9;

    θ_c = 32.781667
    ϕ_c = -79.925

    dst = @. ((θ-θ_c)^2 + (ϕ-ϕ_c)^2 )
    coords = argmin(dst)
    # coords = CartesianIndex(66,50,11)
    ts_ecco_lcl = zeros(Float32,12*length(years))
    ts_ecco_glb = zeros(Float32,12*length(years))
    for yr in years
        for mnth in months
            # println(string(yr)*" "*string(mnth))
            idx = 12*(yr-years[1])+mnth
            fn = homedir()*"/Data/ECCO/v4r4/SSHDYN/"*string(yr)*"/SSHDYN_"*string(yr)*"_"*lpad(string(mnth),2,"0")*".nc"
            sshdyn = ncread(fn,"SSHDYN")[:,:,:,1]
            ts_ecco_glb[idx] = sum(@. area*sshdyn)/sum(area)
            ts_ecco_lcl[idx] = sshdyn[coords]
        end
    end

    ts_ecco_noseas = separate_mss(t,tval,ts_ecco_lcl.-ts_ecco_glb);
    ts_ecco_noseas .-= mean(ts_ecco_noseas);
    return tval,ts_ecco_noseas
end

function movav(t,ts,width)
    @assert isodd(width)
    wrange=convert(Int,width/2-0.5)
    ts_flt = similar(ts)
    @. ts_flt[:] = NaN
    for tstep in 1:length(t)
        tstart=max(1,tstep-wrange)
        tstop=min(length(t),tstep+wrange)
        ts_flt[tstep] = mean(ts[tstart:tstop])
    end
    return ts_flt
end
