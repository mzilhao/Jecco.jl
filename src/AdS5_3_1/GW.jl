#This script will set up and solve the equations, in fourier space, for the gravitational wave
#perturbation h as a function of time. We solve (Dtt+k^2)h(t,k)=T(t,k).

#TODO: The integrator does not evolve in time, always spits dh = 0, fix it.


using FFTW
using LinearAlgebra
using Interpolations
using OrdinaryDiffEq
#Transform of only real functions. Gives back a matrix of 2 component vector for
# the momenta k and the matrix for the mode coefficients fk.
abstract type Parameters end

struct param{T<:Interpolations.GriddedInterpolation, TP<:Real} <: Parameters# where {T<:Interpolations.GriddedInterpolation, TP<:Real}
    f  :: T
    x  :: Array{TP,1}
    y  :: Array{TP,1}
    kx :: Array{TP,2}
    ky :: Array{TP,2}
    dt :: TP
end

function get_evol_variable(u::Array{T,2}, u_t::Array{T,2}) where {T<:Complex}
    N      = length(u)
    u_evol = im.*zeros(N)
    u_evol = reshape(u, N)
    append!(u_evol, reshape(u_t, N))

    u_evol
end

function get_tensors(u_evol::Array{T,1}, Nkx::N, Nky::N) where {T<:Complex, N<:Int}
    Ntot     = length(u_evol)
    idxsplit = Int(Ntot/2)
    h        = reshape(u_evol[1:idxsplit],Nkx,Nky)
    h_t      = reshape(u_evol[idxsplit+1:end],Nkx,Nky)

    h, h_t
end

function Fourier_Transform_2D(f::Array{T,2}) where {T<:Real}
    Nx   = length(f[:,1])
    Ny   = length(f[1,:])
    plan = plan_rfft(f)

    fk   = 1/(Nx*Ny).*(plan * f)
end


function initial_conditions(ff::VEVTimeSeries, kx::Array{T,2}, ky::Array{T,2}) where {T<:Real}
    Tk  = Fourier_Transform_2D(ff[1,:,:])
    Nkx = length(kx[:,1])
    Nky = length(ky[1,:])

    tol = 10^-10
    indeces = findall(abs.(Tk).<tol)
    for i in indeces
        Tk[i] = 0.0
    end


#We can paralelize this using for loops and Threads.@threads macro
#Probably is going to be faster
    h_t    = im.*zeros(Nkx,Nky)
    h      = Tk./(kx.^2 + ky.^2)
    h[1,1] = Tk[1,1]

    get_evol_variable(h,h_t)
end


function rhs(h_evol::Array{T,1}, param::Parameters, t::TP) where {T<:Complex, TP<:Real}
    f_inter = param.f(t,param.x,param.y)
    kx      = param.kx
    ky      = param.ky
    Nkx     = length(kx[:,1])
    Nky     = length(ky[1,:])
    Tk      = Fourier_Transform_2D(f_inter)
    tol     = 10^-10
    indeces = findall(abs.(Tk).<tol)
    for i in indeces
        Tk[i] = 0.0
    end
    h, h_t    = get_tensors(h_evol,Nkx,Nky)
    dh        = h_t
    dh_t      = Tk - (kx.^2 + ky.^2) .*h
    dh_t[1,1] = 0.0

    get_evol_variable(dh,dh_t)
end

#Computes the component associated to field in the bounbdary stress tensor.
#E.g. field=energy gives us h00 as solution.
function solve_GW(dirname::String, field::Symbol; dt::T, alg, option::String="static") where {T<:Real}
    ff       = VEVTimeSeries(dirname, field)
    tt, x, y = get_coords(ff,:,:,:)
    tspan    = (tt[1],tt[end])
    dx       = x[2]-x[1]
    dy       = y[2]-y[1]
    kxx      = 2π.*rfftfreq(length(x), 1/dx)
    kyy      = 2π.*fftfreq(length(y), 1/dy)
    Nkx      = length(kxx)
    Nky      = length(kyy)
    kx       = zeros(Nkx,Nky)
    ky       = zeros(Nkx,Nky)
    for j in eachindex(kyy)
        for i in eachindex(kxx)
            kx[i,j] = kxx[i]
            ky[i,j] = kyy[j]
        end
    end
    h0_evol  = initial_conditions(ff,kx,ky)
    ff_inter = interpolate((tt,x,y,), ff[:,:,:], Gridded(Linear()))
    prm      = param{typeof(ff_inter),typeof(x[1])}(ff_inter,x,y,kx,ky,dt)

    problem    = ODEProblem(rhs, h0_evol, tspan, prm)
    integrator = init(problem, alg, save_everystep=false, dt=dt, adaptive=false)

    h   = im.*zeros(Int(floor((tspan[2]-tspan[1])/dt+1)),Nkx,Nky)
    h_t = similar(h)
    Tk  = similar(h)
    n = 0
    for (u, t) in tuples(integrator)
        n += 1
        h[n,:,:], h_t[n,:,:] = get_tensors(u,Nkx,Nky)
        Tk[n,:,:]            = Fourier_Transform_2D(ff_inter(t,x,y))
        umax  = maximum(abs.(u[2:Nkx*Nky]))
        dumax = maximum(abs.(u[Nkx*Nky+2:end]))

        println("t = $t")
        println("\n")
        println("max |h| = $umax")
        println("max |h_t| = $dumax")
        println("------------------------------------------------------------")
        println("\n")
    end
    h, h_t, Tk, kxx, kyy
end
