# ---------------------------------------------------------------------------------
# Read adjoint sensitivities and interpolate onto regular lat/lon grid for plotting
# Also save sum of sensitivities squared for each lag
# ---------------------------------------------------------------------------------
using Plots
using NCDatasets
using LinearAlgebra
using NetCDF
using Dates
using Statistics
using PyCall
using DelimitedFiles

function main()
    # Read all sensitivities
    processes = ["empmr","qnet","tauu","tauv"]
    ϕ = [0.25:0.5:359.75...]
    θ = [-89.75:0.5:89.75...]

    tval = [2003 + 1/12 - 1/24 : 1/12 : 2004 + 1/12...]
    mnth_array = zeros(Int,length(tval),2)
    mnth_array[:,1] = floor.(tval)
    mnth_array[:,2] = round.(rem.(tval,1) .* 12 .+0.5)

    # time_orig = zeros(Date,107)
    # proc_tot = zeros(107,4)

    fn_save = homedir()*"/Data/ECCO/RISE/results/SanDiego/sens_jan_interp.nc"
    fh = Dataset(fn_save,"c")
    defDim(fh,"time",length(tval))
    defDim(fh,"lon",size(ϕ,1))
    defDim(fh,"lat",size(θ,1))
    defVar(fh,"lon",ϕ,("lon",),fillvalue=-9999,deflatelevel=5)
    defVar(fh,"lat",θ,("lat",),fillvalue=-9999,deflatelevel=5)
    defVar(fh,"time",tval,("time",),fillvalue=-9999,deflatelevel=5)

    sq2 = x -> x^2 
    for (proc_num,proc_name) in enumerate(processes)
        sensitivity = read_sensitivity(proc_name)
        # Total sensitivity
        # proc_tot[:,proc_num] = mapreduce(sq2,+,sensitivity[1].sens,dims=(1,2,3))[1,1,1,:]
        # time_orig[:] = sensitivity[1].time
        ϕ,θ,field_interp = interpolate_2d(sensitivity[1].sens)
        field_interp_mnth = zeros(Float32,length(ϕ),length(θ),length(tval));


        month_sens = [Month(sensitivity[1].time[i]).value for i in 1:length(sensitivity[1].time)]
        year_sens = [Year(sensitivity[1].time[i]).value for i in 1:length(sensitivity[1].time)]
        for mnth in 1:size(mnth_array,1)
            acc_idx = @. (year_sens == mnth_array[mnth,1]) & (month_sens == mnth_array[mnth,2])
            field_interp_mnth[:,:,mnth] = mean(field_interp[:,:,acc_idx],dims=3)[:,:,1]
        end

        field_interp_mnth .*= 1000

        defVar(fh,proc_name,-field_interp_mnth,("lon","lat","time"),deflatelevel=2)
    end
    close(fh)

    # Save total process variance 
    # proc_norm = proc_tot./maximum(proc_tot,dims=1)
    # t_lag = [(time_orig[i] .- Date(2004,6,1)).value for i in 1:length(time_orig)]
    # for (proc_num,proc_name) in enumerate(processes)
    #     fn_save = homedir()*"/Projects/2020_ECCO_adjoint/GMT/total_sens/"*proc_name*".txt"
    #     arr = [t_lag proc_norm[:,proc_num]]
    #     writedlm(fn_save,arr,';')
    # end
    return
end

function read_sensitivity(proc_name)
    println("    Reading adjoint sensitivity...")
    sensitivity = Array{Sensitivity,1}(undef,12)
    dir_sensitivities = homedir()*"/Data/ECCO/RISE/adjoint_sensitivities/SanDiego/"
    yr=2004
    Threads.@threads for mnth in 1:12
        filename = dir_sensitivities*lpad(string(mnth),2,"0")*lpad(string(yr-2000),2,"0")*"/adxx_"*proc_name*".0000000129_dim.data";
        sens_lcl = read_sensitivity_file(filename)
        ts_final = sens_ts_final(mnth)
        sens_ts = ts_final .- [size(sens_lcl,4)-1:-1:0...] .* Week(1) .- Day(3)
        sens_smooth  = similar(sens_lcl)
        sens_smooth[:,:,:,1] = sens_lcl[:,:,:,1]
        for t=2:size(sens_smooth,4)
            @. sens_smooth[:,:,:,t] = 0.5*(sens_lcl[:,:,:,t-1] + sens_lcl[:,:,:,t])
        end
        sensitivity[mnth] = Sensitivity(sens_ts,sens_smooth)
    end
    return sensitivity
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

struct Sensitivity
    time::Array{Date,1}
    sens::Array{Float32,4}
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

    field_interp = zeros(size(ϕ,1),size(θ,1),size(field_ecco,4));
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
    return(ϕ,θ,field_interp)
end


