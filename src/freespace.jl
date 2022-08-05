struct FreeSpaceGreensFunction2D{T} <: AbstractGreensFunction2D
    wavenumber::T
end
const FreeSpace = FreeSpaceGreensFunction2D

function evaluate(G::FreeSpace, x::SVec2, y::SVec2)
    k, r = wavenumber(G), norm(x - y)
    return -(im / 4) * hankelh1(0, k * r)
end

function evaluate(G::FreeSpace, order::SymOrInt, x::SVec2, y::SVec2)
    k, r = wavenumber(G), norm(x - y)

    @assert order in (0, :log, 1, 2)

    if order == 0
        if x == y
            return (log(k / 2) + γ) / 2 / π - im / 4
        else
            return G(x,y) - G(:log,x,y)*cauchyweight2(:log,x,y)
        end
    end
    order == :log && return besselj0(k * r) / 2
    order == 1 && return zero(r)
    order == 2 && return zero(r)
end

function evaluate(G::FreeSpace, x::SVec2, y::SVec2, dx::SVec2)
    k, r = wavenumber(G), norm(x - y)

    return (im * k / 4) * hankelh1_1(k * r) * cosangle(x - y, dx)
end

function evaluate(G::FreeSpace, order::SymOrInt, x::SVec2, y::SVec2, dx::SVec2)
    k, r = wavenumber(G), norm(x - y)

    return (im * k / 4) * hankelh1_1(k * r, order) * cosangle(x - y, dx)
end

function evaluate(G::FreeSpace, x::SVec2, y::SVec2, dx::SVec2, dy::SVec2)
    k, r = wavenumber(G), norm(x - y)

    u = cosangle(x - y, dx) * cosangle(x - y, dy)
    v = sum(dx .* dy)

    return (im * k^2 / 4) * (u * hankelh1_2(k * r) - v * hinc(k * r))
end

function evaluate(G::FreeSpace, order::SymOrInt, x::SVec2, y::SVec2, dx::SVec2, dy::SVec2)
    k, r = wavenumber(G), norm(x - y)

    u = cosangle(x - y, dx) * cosangle(x - y, dy)
    v = sum(dx .* dy)

    return (im * k^2 / 4) * (u * hankelh1_2(k * r, order) - v * hinc(k * r, order))
end
