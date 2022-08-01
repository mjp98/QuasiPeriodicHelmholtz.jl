abstract type AbstractWaveguideGreens2D <: AbstractGreensFunction2D end

struct WaveguideNeumann{T,S<:Number} <: AbstractWaveguideGreens2D
    G::T
    d::S
end
struct WaveguideDirichlet{T,S<:Number} <: AbstractWaveguideGreens2D
    G::T
    d::S
end
# struct WaveguideMixed{T,S<:Number} <: AbstractWaveguideGreens2D
#     G::T
#     d::S
# end

greens(G::AbstractWaveguideGreens2D) = G.G

function evaluate(G::WaveguideNeumann,x::SVec2,ξ::SVec2;M::Int=256)
    g = greens(G)
    flipx = Diagonal(@SVector [-1,1])
    return evaluate(g,x,ξ;M) + evaluate(g,x,flipx*ξ;M)
end
function evaluate(G::WaveguideDirichlet,x::SVec2,ξ::SVec2;M::Int=256)
    g = greens(G)
    flipx = Diagonal(@SVector [-1,1])
    return evaluate(g,x,ξ;M) - evaluate(g,x,flipx*ξ;M)
end

# function evaluate(G::WaveguideMixed,x::SVec2,ξ::SVec2;M::Int=256)
#     g = Greens(G)
#     flipx = Diagonal(@SVector [-1,1])
#     return evaluate(g,x,ξ;M) - evaluate(g,x,flipx*ξ;M)
# end
