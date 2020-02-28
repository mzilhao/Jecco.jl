
abstract type PotentialType end
abstract type Const  <: PotentialType end
abstract type Square <: PotentialType end

struct Potential{fType,fpType,T<:PotentialType} <: Function
    f  :: fType
    fp :: fpType
end
Potential{T}(f,fp) where {T} = Potential{typeof(f),typeof(fp),T}(f,fp)

(ff::Potential)(phi)  = ff.f(phi)
∂(ff::Potential)      = ff.fp

Potential{Const}()  = Potential{Const}(phi -> phi, phi-> 1.0)
Potential{Square}() = Potential{Square}(phi -> -1.0 + 0.5 * phi*phi, phi -> phi)

function Potential(par_base::ParamBase)
    if par_base.which_potential == "square"
        Potential{Square}()
    elseif par_base.which_potential == "const"
        Potential{Const}()
    else
        error("Unknown potential.")
    end
end
