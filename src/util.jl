const RealOrDualReal = Union{Real,Dual{<:Real}}
const SVec2 = SVector{2}
const SymOrInt = Union{Symbol,Integer}

logabs(z) = log(norm(z))

function cauchyweight2(order ,x)
    order == :log && return logabs(x) / π
    order == 0 && return 1
    return abs(x)^(-order) / π
end
cauchyweight2(order, x, y) = cauchyweight2(order,y-x)

derivative(f, z) = epsilon(f(Dual(z, 1)))

function cosangle(x::SVec2, dx::SVec2)
    return iszero(x) ? zero(eltype(x)) : dot(x, dx) / norm(x)
end
