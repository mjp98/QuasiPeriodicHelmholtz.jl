
abstract type AbstractGreensFunction2D <: Function end
const Greens = AbstractGreensFunction2D


wavenumber(G::Greens) = G.wavenumber
(G::Greens)(args...;kwargs...) = evaluate(G, args...;kwargs...)

# function evaluate(G::Greens, x::Complex, args...;kwargs...)
#     return evaluate(G,SVector(reim(x)...),(SVector(reim(z)...) for z in args...)...;kwargs...)
# end




function evaluate(G::Greens, x::Complex, y::Complex;kwargs...)
    evaluate(G,SVector(reim(x)...),SVector(reim(y)...);kwargs...)
end

function evaluate(G::Greens, x::Complex, y::Complex, dx::Complex;kwargs...)
    evaluate(G,SVector(reim(x)...),SVector(reim(y)...),SVector(reim(dx)...);kwargs...)
end

function evaluate(G::Greens, x::Complex, y::Complex, dx::Complex,dy::Complex;kwargs...)
    evaluate(G,SVector(reim(x)...),SVector(reim(y)...),SVector(reim(dx)...),SVector(reim(dy)...);kwargs...)
end

function evaluate(G::Greens, order::SymOrInt,x::Complex, y::Complex;kwargs...)
    evaluate(G,order,SVector(reim(x)...),SVector(reim(y)...);kwargs...)
end

function evaluate(G::Greens, order::SymOrInt,x::Complex, y::Complex, dx::Complex;kwargs...)
    evaluate(G,order,SVector(reim(x)...),SVector(reim(y)...),SVector(reim(dx)...);kwargs...)
end

function evaluate(G::Greens, order::SymOrInt,x::Complex, y::Complex, dx::Complex,dy::Complex;kwargs...)
    evaluate(G,order,SVector(reim(x)...),SVector(reim(y)...),SVector(reim(dx)...),SVector(reim(dy)...);kwargs...)
end


# function evaluate(G::Greens, order::SymOrInt, x::Complex, args...;kwargs...)
#     return evaluate(G,order, SVector(reim(x)...),(SVector(reim(z)...) for z in args...)...;kwargs...)
# end
