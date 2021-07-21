#This script will set up and solve the equations, in fourier space, for the gravitational wave
#perturbation h as a function of time. We solve (Dtt+k^2)h(t,k)=T(t,k).

#TODO: The integrator does not evolve in time, always spits dh = 0, fix it.

import Base.Threads.@threads
using FFTW
using LinearAlgebra
using Interpolations
using OrdinaryDiffEq
#Transform of only real functions. Gives back a matrix of 2 component vector for
# the momenta k and the matrix for the mode coefficients fk.
abstract type Parameters end

struct param{T<:Interpolations.GriddedInterpolation, TP<:Real} <: Parameters
    px  :: T
    pxy :: T
    py  :: T
    pz  :: T
    x   :: Array{TP,1}
    y   :: Array{TP,1}
    kx  :: Array{TP,2}
    ky  :: Array{TP,2}
    dt  :: TP
    tol :: TP
end

function get_evol_variable(u::Array{T,3}, u_t::Array{T,3}) where {T<:Complex}
    N      = length(u)
    u_evol = im.*zeros(N)
    u_evol = reshape(u, N)
    append!(u_evol, reshape(u_t, N))

    u_evol
end

function get_tensors(u_evol::Array{T,1}, Nkx::N, Nky::N) where {T<:Complex, N<:Int}
    Ntot     = length(u_evol)
    idxsplit = Int(Ntot/2)
    h        = reshape(u_evol[1:idxsplit],Nkx,Nky,4)
    h_t      = reshape(u_evol[idxsplit+1:end],Nkx,Nky,4)

    h, h_t
end

function Fourier_Transform_2D(f::Array{T,2}) where {T<:Real}
    Nx   = length(f[:,1])
    Ny   = length(f[1,:])
    plan = plan_rfft(f)

    fk   = 1/(Nx*Ny).*(plan * f)
end

function initial_conditions(ff::VEVTimeSeries, kx::Array{T,2}, ky::Array{T,2}) where {T<:Real}
    #Tk  = Fourier_Transform_2D(ff[1,:,:])
    Nkx = length(kx[:,1])
    Nky = length(ky[1,:])


    h_t    = im.*zeros(Nkx, Nky, 4)
    h      = im.*zeros(Nkx, Nky, 4)


    get_evol_variable(h,h_t)
end


function rhs(h_evol::Array{T,1}, param::Parameters, t::TP) where {T<:Complex, TP<:Real}
    tol       = param.tol
    px_inter  = param.px(t,param.x,param.y)
    pxy_inter = param.pxy(t,param.x,param.y)
    py_inter  = param.py(t,param.x,param.y)
    pz_inter  = param.pz(t,param.x,param.y)
    kx        = param.kx
    ky        = param.ky
    Nkx       = length(kx[:,1])
    Nky       = length(ky[1,:])
    px        = Fourier_Transform_2D(px_inter)
    pxy       = Fourier_Transform_2D(pxy_inter)
    py        = Fourier_Transform_2D(py_inter)
    pz        = Fourier_Transform_2D(pz_inter)

    h, h_t    = get_tensors(h_evol, Nkx, Nky)
    dh        = similar(h)
    dh_t      = similar(h_t)

#TODO: Do proper paralesization, @threads does not work like this as you are using
#the same name for the variable kkx, kky... and different threads will mess up
#with it.
#Also it might be better to solve all kx and ky as matrix equation at once
#and do a loop over the 4 components that we have to solve.
#Bear in mind that we will get to runs with many points in x and y.
    @time @fastmath @inbounds @threads for j in 1:Nky
        for i in 1:Nkx
            kkx  = kx[i,j]
            kky  = ky[i,j]
            kkx2 = kkx^2
            kky2 = kky^2
            k2   = kkx^2+kky^2
            ppx  = px[i,j]
            ppxy = pxy[i,j]
            ppy  = py[i,j]
            ppz  = pz[i,j]

            M = 0.5/k2*(kky2*ppx+kkx2*ppy+2*kkx*kky*ppxy-k2*ppz).*[kky^2/k2, kkx*kky/k2, kkx^2/k2, -1]
            indices = findall(abs.(M) .< tol)
            for i in indices
                M[i] = 0.0
            end
            dh[i,j,:]   = h_t[i,j,:]
            dh_t[i,j,:] = -k2.*h[i,j,:]-M
        end
    end
    ppx         = px[1,1]
    ppxy        = pxy[1,1]
    ppy         = py[1,1]
    ppz         = pz[1,1]
    trT         = ppx + ppy + ppz
    M           = [ppx, ppxy, ppy, ppz] -1/3*trT.*[1,0,1,1]
    indices     = findall(abs.(M) .< tol)
    for i in indices
        M[i]    = 0.0
    end
    dh_t[1,1,:] = -M

    get_evol_variable(dh,dh_t)
end

#Computes the component associated to field in the bounbdary stress tensor.
#E.g. field=energy gives us h00 as solution.
function solve_GW(dirname::String; dt::T = 0.0, alg, option::String="static", tol::T=10^-10) where {T<:Real}
    px       = VEVTimeSeries(dirname, :px)
    pxy      = VEVTimeSeries(dirname, :pxy)
    py       = VEVTimeSeries(dirname, :py)
    pz       = VEVTimeSeries(dirname, :pz)
    tt, x, y = get_coords(px,:,:,:)
    tspan    = (tt[1],tt[end])
    if dt == 0.0 dt = tt[2]-tt[1] end
    dx       = x[2]-x[1]
    dy       = y[2]-y[1]
    kxx      = 2π.*rfftfreq(length(x), 1/dx)
    kyy      = 2π.*fftfreq(length(y), 1/dy)
    Nkx      = length(kxx)
    Nky      = length(kyy)
    kx       = zeros(Nkx,Nky)
    ky       = zeros(Nkx,Nky)

    @time @fastmath @inbounds @threads for j in eachindex(kyy)
        for i in eachindex(kxx)
            kx[i,j] = kxx[i]
            ky[i,j] = kyy[j]
        end
    end

    h0_evol   = initial_conditions(px,kx,ky)
    px_inter  = interpolate((tt,x,y,), px[:,:,:], Gridded(Linear()))
    pxy_inter = interpolate((tt,x,y,), pxy[:,:,:], Gridded(Linear()))
    py_inter  = interpolate((tt,x,y,), py[:,:,:], Gridded(Linear()))
    pz_inter  = interpolate((tt,x,y,), pz[:,:,:], Gridded(Linear()))
    ff_inter  = [px_inter, pxy_inter, py_inter, pz_inter]
    prm       = param{typeof(px_inter),typeof(x[1])}(px_inter,pxy_inter,py_inter,pz_inter,
                                                            x,y,kx,ky,dt,tol)
    problem    = ODEProblem(rhs, h0_evol, tspan, prm)
    integrator = init(problem, alg, save_everystep=false, dt=dt, adaptive=false)

    h   = im.*zeros(Int(floor((tspan[2]-tspan[1])/dt+1)),Nkx,Nky,4)
    h_t = similar(h)
    t_2 = zeros(length(h[:,1,1,1]))
    pxk = im.*zeros(length(t_2), Nkx, Nky)
    pyk = im.*zeros(length(t_2), Nkx, Nky)
    pzk = im.*zeros(length(t_2), Nkx, Nky)
    uu = im.*zeros(length(t_2),length(h0_evol))
    n = 0
    for (u, t) in tuples(integrator)
        n += 1
        t_2[n]                   = t
        h[n,:,:,:], h_t[n,:,:,:] = get_tensors(u,Nkx,Nky)
        pxk[n,:,:]   = Fourier_Transform_2D(prm.px(t,prm.x,prm.y))
        pyk[n,:,:]   = Fourier_Transform_2D(prm.py(t,prm.x,prm.y))
        pzk[n,:,:]   = Fourier_Transform_2D(prm.pz(t,prm.x,prm.y))
        umax  = maximum(abs.(u[2:Nkx*Nky]))
        dumax = maximum(abs.(u[Nkx*Nky+2:end]))
        indices = findfirst(abs.(h_t[n,:,:,:]) .== dumax)

        println("t = $t")
        println("\n")
        println("max |h| = $umax")
        println("max |h_t| = $dumax")
        println("index = $indices")
        println("------------------------------------------------------------")
        println("\n")
    end
    h, h_t, t_2, kxx, kyy, pxk, pyk, pzk
end
