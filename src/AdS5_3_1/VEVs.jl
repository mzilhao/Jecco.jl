
"""
    compute_energy(a4, phi2, phi0, phiM2)

Returns the normalized energy density ``ϵ / ϕ₀⁴``.
"""
function compute_energy(a4, phi2, phi0, oophiM2)
    if phi0 < 1.e-9
        return -3//4 * a4
    else
        return -3//4 * a4 / phi0^4 - phi2 / phi0^3 .- 1//4 * oophiM2 .+ 7//36
    end
end

"""
    compute_flux(f2)

Returns the normalized flux ``J / ϕ₀⁴``.
"""
function compute_flux(f2, phi0)
    if phi0 < 1.e-9
        return f2
    else
        return f2 / phi0^4
    end
end

"""
    compute_pressure_x(a4, b14, b24, phi2, phi0, oophiM2)

Returns the normalized pressure in the x direction ``pₓ / ϕ₀⁴``.
"""
function compute_pressure_x(a4, b14, b24, phi2, phi0, oophiM2)
    if phi0 < 1.e-9
        return -a4/4 - b14 - b24
    else
        return -a4/(4*phi0^4) - b14/phi0^4 - b24/phi0^4 + phi2/(3*phi0^3) .+
            (-5//108 + 1//4 * oophiM2)
    end
end

"""
    compute_pressure_y(a4, b14, b24, phi2, phi0, oophiM2)

Returns the normalized pressure in the y direction ``p_y / ϕ₀⁴``.
"""
function compute_pressure_y(a4, b14, b24, phi2, phi0, oophiM2)
    if phi0 < 1.e-9
        return -a4/4 + b14 - b24
    else
        return -a4/(4*phi0^4) + b14/phi0^4 - b24/phi0^4 + phi2/(3*phi0^3) .+
            (-5//108 + 1//4*oophiM2)
    end
end

"""
    compute_pressure_z(a4, b24, phi2, phi0, oophiM2)

Returns the normalized pressure in the z direction ``p_z / ϕ₀⁴``.
"""
function compute_pressure_z(a4, b24, phi2, phi0, oophiM2)
    if phi0 < 1.e-9
        return -a4/4 + 2*b24
    else
        return -a4/(4*phi0^4) + 2*b24/phi0^4 + phi2/(3*phi0^3) .+
            (-5//108 + 1//4*oophiM2)
    end
end

"""
    compute_pressure_xy(a4, b24, phi2, phi0, oophiM2)

Returns the normalized pressure in the xy direction ``p_xy / ϕ₀⁴``.
"""
function compute_pressure_xy(g4, phi0)
    if phi0 < 1.e-9
        return -g4
    else
        return -g4 / phi0^4
    end
end

"""
    compute_Ophi(phi2, phi0, oophiM2)

Returns the normalized scalar operator VEV ``O_ϕ / ϕ₀³``.
"""
function compute_Ophi(phi2, phi0, oophiM2)
    if phi0 < 1.e-9
        return -2 * phi2
    else
        return -2 * phi2 / phi0^3 .+ (1//3 - oophiM2)
    end
end


function get_energy(ts::OpenPMDTimeSeries; it::Int, verbose::Bool=false)
    a4,   chart = get_field(ts, it=it, field="a4", verbose=verbose)
    phi2, chart = get_field(ts, it=it, field="phi2", verbose=verbose)

    phi0 = try
        ts.params["phi0"]
    catch e
        if isa(e, KeyError)
            0.0   # if "phi0" is not found in the params Dict, set phi0 = 0
        else
            throw(e)
        end
    end

    oophiM2 = try
        ts.params["oophiM2"]
    catch e
        if isa(e, KeyError)
            0.0   # if "oophiM2" is not found in the params Dict, set oophiM2 = 0
        else
            throw(e)
        end
    end

    en = compute_energy(a4, phi2, phi0, oophiM2)

    en, chart
end

function get_Jx(ts::OpenPMDTimeSeries; it::Int, verbose::Bool=false)
    fx2, chart = get_field(ts, it=it, field="fx2", verbose=verbose)

    phi0 = try
        ts.params["phi0"]
    catch e
        if isa(e, KeyError)
            0.0   # if "phi0" is not found in the params Dict, set phi0 = 0
        else
            throw(e)
        end
    end

    Jx = compute_flux(fx2, phi0)

    Jx, chart
end

function get_Jy(ts::OpenPMDTimeSeries; it::Int, verbose::Bool=false)
    fy2, chart = get_field(ts, it=it, field="fy2", verbose=verbose)

    phi0 = try
        ts.params["phi0"]
    catch e
        if isa(e, KeyError)
            0.0   # if "phi0" is not found in the params Dict, set phi0 = 0
        else
            throw(e)
        end
    end

    Jy = compute_flux(fy2, phi0)

    Jy, chart
end

function get_px(ts::OpenPMDTimeSeries; it::Int, verbose::Bool=false)
    a4,   chart = get_field(ts, it=it, field="a4", verbose=verbose)
    b14,  chart = get_field(ts, it=it, field="b14", verbose=verbose)
    b24,  chart = get_field(ts, it=it, field="b24", verbose=verbose)
    phi2, chart = get_field(ts, it=it, field="phi2", verbose=verbose)

    phi0 = try
        ts.params["phi0"]
    catch e
        if isa(e, KeyError)
            0.0   # if "phi0" is not found in the params Dict, set phi0 = 0
        else
            throw(e)
        end
    end

    oophiM2 = try
        ts.params["oophiM2"]
    catch e
        if isa(e, KeyError)
            0.0   # if "oophiM2" is not found in the params Dict, set oophiM2 = 0
        else
            throw(e)
        end
    end

    px = compute_pressure_x(a4, b14, b24, phi2, phi0, oophiM2)

    px, chart
end

function get_py(ts::OpenPMDTimeSeries; it::Int, verbose::Bool=false)
    a4,   chart = get_field(ts, it=it, field="a4", verbose=verbose)
    b14,  chart = get_field(ts, it=it, field="b14", verbose=verbose)
    b24,  chart = get_field(ts, it=it, field="b24", verbose=verbose)
    phi2, chart = get_field(ts, it=it, field="phi2", verbose=verbose)

    phi0 = try
        ts.params["phi0"]
    catch e
        if isa(e, KeyError)
            0.0   # if "phi0" is not found in the params Dict, set phi0 = 0
        else
            throw(e)
        end
    end

    oophiM2 = try
        ts.params["oophiM2"]
    catch e
        if isa(e, KeyError)
            0.0   # if "oophiM2" is not found in the params Dict, set oophiM2 = 0
        else
            throw(e)
        end
    end

    py = compute_pressure_y(a4, b14, b24, phi2, phi0, oophiM2)

    py, chart
end

function get_pz(ts::OpenPMDTimeSeries; it::Int, verbose::Bool=false)
    a4,   chart = get_field(ts, it=it, field="a4", verbose=verbose)
    b24,  chart = get_field(ts, it=it, field="b24", verbose=verbose)
    phi2, chart = get_field(ts, it=it, field="phi2", verbose=verbose)

    phi0 = try
        ts.params["phi0"]
    catch e
        if isa(e, KeyError)
            0.0   # if "phi0" is not found in the params Dict, set phi0 = 0
        else
            throw(e)
        end
    end

    oophiM2 = try
        ts.params["oophiM2"]
    catch e
        if isa(e, KeyError)
            0.0   # if "oophiM2" is not found in the params Dict, set oophiM2 = 0
        else
            throw(e)
        end
    end

    pz = compute_pressure_z(a4, b24, phi2, phi0, oophiM2)

    pz, chart
end

function get_pxy(ts::OpenPMDTimeSeries; it::Int, verbose::Bool=false)
    g4,   chart = get_field(ts, it=it, field="g4", verbose=verbose)

    phi0 = try
        ts.params["phi0"]
    catch e
        if isa(e, KeyError)
            0.0   # if "phi0" is not found in the params Dict, set phi0 = 0
        else
            throw(e)
        end
    end

    pxy = compute_pressure_xy(g4, phi0)

    pxy, chart
end

function get_Ophi(ts::OpenPMDTimeSeries; it::Int, verbose::Bool=false)
    phi2,  chart = get_field(ts, it=it, field="phi2", verbose=verbose)

    phi0 = try
        ts.params["phi0"]
    catch e
        if isa(e, KeyError)
            0.0   # if "phi0" is not found in the params Dict, set phi0 = 0
        else
            throw(e)
        end
    end

    oophiM2 = try
        ts.params["oophiM2"]
    catch e
        if isa(e, KeyError)
            0.0   # if "oophiM2" is not found in the params Dict, set oophiM2 = 0
        else
            throw(e)
        end
    end

    Ophi = compute_Ophi(phi2, phi0, oophiM2)

    Ophi, chart
end

#=
Computes the fluid velocity and the local frame (diagonal) stress tensor. p3 is just pz. I think that
ut will always be 1, if not we can always hand ux/ut and uy/ut, which diretly are the speeds measured
in the LAB frame, i.e. dx/dt and dy/dt.
=#

function compute_local_VEVs(e, Jx, Jy, px, py, pxy)
    T = [e -Jx -Jy; -Jx px pxy; -Jy pxy py]
    values, vectors = eigen(T)
    norms = -vectors[1,:].^2 + vectors[2,:].^2 + vectors[3,:].^2
    index = findfirst(==(1), norms.<0)
    ut = sign(vectors[1,index])*vectors[1,index]/(-norms[index])^0.5
    ux = sign(vectors[1,index])*vectors[2,index]/(-norms[index])^0.5
    uy = sign(vectors[1,index])*vectors[3,index]/(-norms[index])^0.5
    e_local = values[index]
    p1_local = values[1:end .!=index][1]
    p2_local = values[1:end .!=index][2]

    return ut, ux, uy, e_local, p1_local, p2_local
end
