using HDF5

# TODO: make phi0 and oophiM2 part of the hdf5 attributes?
function convert_to_mathematica(dirname::String; outfile::String="data_mathematica.h5",
                                phi0, oophiM2)
    ts = OpenPMDTimeSeries(dirname, prefix="boundary_")

    iterations = ts.iterations

    Nt = length(iterations)
    t  = zeros(Nt)

    it = 0
    en0, chart = get_energy(ts, it=it, phi0=phi0, oophiM2=oophiM2)

    _, x, y = chart[:]
    Nx = length(x)
    Ny = length(y)

    en   = zeros(Nt,Nx,Ny)
    Jx   = zeros(Nt,Nx,Ny)
    Jy   = zeros(Nt,Nx,Ny)
    px   = zeros(Nt,Nx,Ny)
    py   = zeros(Nt,Nx,Ny)
    pz   = zeros(Nt,Nx,Ny)
    pxy  = zeros(Nt,Nx,Ny)
    Ophi = zeros(Nt,Nx,Ny)

    for (idx,it) in enumerate(iterations)
        en[idx,:,:]   .= get_energy(ts, it=it, phi0=phi0, oophiM2=oophiM2)[1][1,:,:]
        Jx[idx,:,:]   .= get_Jx(ts, it=it, phi0=phi0, oophiM2=oophiM2)[1][1,:,:]
        Jy[idx,:,:]   .= get_Jy(ts, it=it, phi0=phi0, oophiM2=oophiM2)[1][1,:,:]
        px[idx,:,:]   .= get_px(ts, it=it, phi0=phi0, oophiM2=oophiM2)[1][1,:,:]
        py[idx,:,:]   .= get_py(ts, it=it, phi0=phi0, oophiM2=oophiM2)[1][1,:,:]
        pz[idx,:,:]   .= get_pz(ts, it=it, phi0=phi0, oophiM2=oophiM2)[1][1,:,:]
        pxy[idx,:,:]  .= get_pxy(ts, it=it, phi0=phi0, oophiM2=oophiM2)[1][1,:,:]
        Ophi[idx,:,:] .= get_Ophi(ts, it=it, phi0=phi0, oophiM2=oophiM2)[1][1,:,:]

        t[idx] = ts.current_t
    end

    # store in an array suitable for Mathematica
    T_m = zeros(11, Nt*Nx*Ny)

    n = 1
    for i in 1:Nt
        for j in 1:Nx
            for k in 1:Ny
                T_m[:,n] = [t[i] x[j] y[k] en[i,j,k] Jx[i,j,k] Jy[i,j,k] px[i,j,k] pxy[i,j,k] py[i,j,k] pz[i,j,k] Ophi[i,j,k]]
                n += 1
            end
        end
    end

    output   = abspath(dirname, outfile)
    fid      = h5open(output, "w")
    group_st = g_create(fid, "data")
    Jecco.write_dataset(group_st, "VEVs", T_m)
    close(fid)
end
