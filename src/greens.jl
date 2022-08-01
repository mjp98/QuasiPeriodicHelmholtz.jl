
abstract type AbstractGreensFunction2D <: Function end
const Greens = AbstractGreensFunction2D


wavenumber(G::Greens) = G.wavenumber
(G::Greens)(x::SVec2, args...;kwargs...) = evaluate(G, x,args...;kwargs...)

function evaluate(G::Greens, x::Complex, args...;kwargs...)
    return evaluate(G,SVector(reim(x)...),(SVector(reim(z)...) for z in args...)...;kwargs...)
end
