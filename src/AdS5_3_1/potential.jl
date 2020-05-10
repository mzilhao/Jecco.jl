
abstract type Potential end

struct ZeroPotential <: Potential end

UU(phi, ::ZeroPotential)  = 0
UUp(phi, ::ZeroPotential) = 0
